library core3d;

import 'dart:html';
import 'dart:web_gl';
import 'dart:typed_data';
import 'dart:math';


/**
 *
 */
class Mesh
{
  Matrix4 transform;      // transform for this mesh
  int numTris;            // num of triangles to render
  List<Mesh> childMeshes; // list of children

  List<double> _vertData;
  List<int> _idxData;
  ImageElement _diffImg;
  ImageElement _normImg;
  ImageElement _specImg;

  Program _prog;          // compiled shader program for this mesh
  Buffer _vertBuffer;     // uploaded vertices data to GPU
  Buffer _idxBuffer;      // uploaded indices data to GPU
  Texture _diffMap;       // diffuse map
  Texture _normMap;       // normals map
  Texture _specMap;       // specular strength map
  Matrix4 _workingT;      // used for rendering
  int _prepOnState=-1;    // current state id shader programs are compiled on

  static RenderingContext _gl;       // reference to the current rendering context
  static List<PointLight> _lightPts; // global point lights
  static int _stateId=0;             // incremented when user changes render parameters
  static Matrix4 _viewT;             // view transform of the whole scene to be rendered
  static Matrix4 _camT;              // camera transform, inverse of viewT
  static Matrix4 _persT;             // perspective transform matrix

  /**
   * creates new mesh from given geometry data
   */
  Mesh([List<double> vertexData=null,List<int> indexData=null])
  {
    // ----- set default global states if not initialized
    if (_viewT==null || _persT==null) Mesh.setCamera(0.0,0.0,10.0,0.0,0.0,0.0);
    if (_lightPts==null) _lightPts = [new PointLight(0.0,0.0,10.0)];

    // ----- set default for this mesh
    transform = new Matrix4();
    childMeshes = new List<Mesh>();

    setGeometry(vertexData,indexData);
  }//endconstr

  /**
   * upload data and compile shader codes
   */
  void _prepForRender()
  {
    if (_vertData!=null && _idxData!=null) setGeometry(_vertData,_idxData);
    if (_diffImg!=null) setTexture(_diffImg);
    if (_normImg!=null) setNormalMap(_normImg);
    if (_specImg!=null) setSpecularMap(_specImg);
    _compileShaderCodes();
  }//endfunction

  /**
   * compiles the shader program specific to the mesh and upload
   */
  void _compileShaderCodes()
  {
    // ----- construct vertex shader source -------------
    String vertSrc = '''
     attribute vec3 aPosn;
     attribute vec2 aUV;
     attribute vec3 aNorm;
     attribute vec3 aTang;
     
     varying vec3 vPosn;
     varying vec2 vUV;
     varying vec3 vNorm;
     varying vec3 vTang;
     
     uniform mat4 uMVMatrix;
     uniform mat4 uPMatrix;
 
     void main(void) {
         vPosn = (uMVMatrix * vec4(aPosn, 1.0)).xyz;
         gl_Position = uPMatrix * vec4(vPosn, 1.0);
         vUV = vec2(aUV.x,-aUV.y);
         vNorm = (uMVMatrix * vec4(aNorm, 0.0)).xyz;
         vTang = (uMVMatrix * vec4(aTang, 0.0)).xyz;
     }
    ''';

    // ----- construct fragment shader source -----------
    String fragSrc ="precision highp float;\n"+
                    "varying vec3 vPosn;\n"+
                    "varying vec2 vUV;\n"+
                    "varying vec3 vNorm;\n"+
                    "varying vec3 vTang;\n"+
                    "uniform vec3 lPoints["+_lightPts.length.toString()+"];\n"+
                    "uniform vec3 lColors["+_lightPts.length.toString()+"];\n";
    if (_diffMap!=null) fragSrc += "uniform sampler2D diffMap;\n";
    if (_normMap!=null) fragSrc += "uniform sampler2D normMap;\n";
    if (_specMap!=null) fragSrc += "uniform sampler2D specMap;\n";

    fragSrc+= "void main(void) {\n"+
              "   vec3 norm = normalize(vNorm);\n"+
              "   vec4 accu = vec4(0.0,0.0,0.0,0.0);\n"+
              "   vec3 color;\n"+       // working var
              "   vec3 lightDir;\n"+    // working var
              "   float f = 0.0;\n\n";  // working var
    if (_normMap!=null) fragSrc +="   color = 2.0*texture2D(normMap, vUV).xyz - vec3(1.0,1.0,1.0);\n"+
                                  "   vec3 tang = normalize(vTang);\n"+
                                  "   vec3 cota = cross(norm,tang);\n"+
                                  "   norm = color.x*tang + color.y*cota + color.z*norm;\n\n";
    if (_diffMap!=null) fragSrc +="   vec4 texColor = texture2D(diffMap, vUV);\n\n";
    else               fragSrc += "   vec4 texColor = vec4(1.0,1.0,1.0,1.0);\n\n";
    for (int i=0; i<_lightPts.length; i++)
    {
      String si = i.toString();
      fragSrc +="   color = lColors["+si+"] * texColor.xyz;\n"+         // color = lightColor*texColor
                "   lightDir = normalize(lPoints["+si+"] - vPosn);\n"+  // light direction from point
                "   f = max(0.0,dot(lightDir,norm));\n"+                // diffuse strength
                "   accu = max(accu,vec4(color*f, texColor.a));\n";     // diffuse shaded color

      fragSrc +="   f = dot(lightDir,normalize(vPosn-2.0*dot(vPosn,norm)*norm));\n"+ // spec strength
                "   f = max(0.0,f*6.0-5.0);\n"+                         // spec strength with hardness
                "   color = f*lColors["+si+"];\n";                      // reflected light intensity
      if (_specMap!=null) fragSrc +="   color = color*texture2D(diffMap, vUV).xyz;\n";
      fragSrc +="   accu = max(accu,vec4(color, texColor.a));\n";       // specular refl color
    }
    fragSrc += "   gl_FragColor = accu;\n}";
    _prog = _uploadShaderProgram(vertSrc,fragSrc);
    if (_prog!=null) _prepOnState = _stateId;  // specify shader code is built based on this state
  }//endfunction

  /**
   * adds given mesh to this render tree
   */
  void addChild(Mesh m)
  {
    childMeshes.add(m);
  }//endfunction

  /**
   * to return as a flattened list the children and grandchildrens of this mesh, self included
   */
  void _flattenTree(Matrix4 T,List<Mesh> L)
  {
    _workingT = T.mult(transform); // working transform of this mesh
    L.add(this);
    for (int i=childMeshes.length-1; i>-1; i--)
      childMeshes[i]._flattenTree(_workingT,L);
  }//endfunction

  /**
   * set vertices and indices of mesh and upload to buffer if context is available
   */
  void setGeometry([List<double> vertexData=null,List<int> indexData=null])
  {
    _vertData = vertexData;
    _idxData = indexData;

    // ----- upload vertex data to buffer
    if (_vertData!=null)
    {
      if (_idxData==null)
      {
        numTris = _vertData.length~/24;
        _idxData = new List<int>();
        for (int i=0; i<numTris*3; i++)  _idxData.add(i);
      }
      else
        numTris = _idxData.length~/3;

      // ----- calc default normals to data if normals are 0,0,0 ----
      for (int i=0; i<numTris; i++)
      {
         int i0=_idxData[i*3+0]*8;
         int i1=_idxData[i*3+1]*8;
         int i2=_idxData[i*3+2]*8;

         if ((_vertData[i0+5]==0 && _vertData[i0+6]==0 && _vertData[i0+7]==0) ||
             (_vertData[i1+5]==0 && _vertData[i1+6]==0 && _vertData[i1+7]==0) ||
             (_vertData[i2+5]==0 && _vertData[i2+6]==0 && _vertData[i2+7]==0))
         {
           // ----- calculate default normals ------------------------
           double vax = _vertData[i0+0];
           double vay = _vertData[i0+1];
           double vaz = _vertData[i0+2];
           double px = _vertData[i1+0] - vax;
           double py = _vertData[i1+1] - vay;
           double pz = _vertData[i1+2] - vaz;
           double qx = _vertData[i2+0] - vax;
           double qy = _vertData[i2+1] - vay;
           double qz = _vertData[i2+2] - vaz;
           // normal by determinant
           double nx = py*qz-pz*qy; //  unit normal x for the triangle
           double ny = pz*qx-px*qz; //  unit normal y for the triangle
           double nz = px*qy-py*qx; //  unit normal z for the triangle
           double nl = sqrt(nx*nx+ny*ny+nz*nz);
           nx/=nl; ny/=nl; nz/=nl;
           _vertData[i0+5]=nx; _vertData[i0+6]=ny; _vertData[i0+7]=nz;
           _vertData[i1+5]=nx; _vertData[i1+6]=ny; _vertData[i1+7]=nz;
           _vertData[i2+5]=nx; _vertData[i2+6]=ny; _vertData[i2+7]=nz;
         }
      }
    }//endif

    if (_vertData==null || _idxData==null || _gl==null)
    {
      numTris=0;
      _prog=null;
      return;
    }

    _idxBuffer = _gl.createBuffer();
    _gl.bindBuffer(ELEMENT_ARRAY_BUFFER, _idxBuffer);
    _gl.bufferDataTyped(ELEMENT_ARRAY_BUFFER,new Uint16List.fromList(_idxData),STATIC_DRAW);

    _vertBuffer = _gl.createBuffer();
    _gl.bindBuffer(ARRAY_BUFFER, _vertBuffer);
    _gl.bufferDataTyped(ARRAY_BUFFER,new Float32List.fromList(calcTangentBasis(_idxData,_vertData)),STATIC_DRAW);
  }//endfunction

  /**
   * sets the diffuse texture for this mesh
   */
  void setTexture(ImageElement img,[bool propagate=false])
  {
    _diffImg = img;
    _diffMap = createMipMapTexture(img);
    _compileShaderCodes();
    if (propagate)
      for (int i=childMeshes.length-1; i>-1; i--)
        childMeshes[i].setTexture(img,propagate);
  }//endfunction

  /**
   * sets the normal map for this mesh
   */
  void setNormalMap(ImageElement img,[bool propagate=false])
  {
    _normImg = img;
    _normMap = createMipMapTexture(img);
    _compileShaderCodes();
    if (propagate)
      for (int i=childMeshes.length-1; i>-1; i--)
        childMeshes[i].setNormalMap(img,propagate);
  }//endfunction

  /**
   * sets the specular strength map for this mesh
   */
  void setSpecularMap(ImageElement img,[bool propagate=false])
  {
    _specImg = img;
    _specMap = createMipMapTexture(img);
    _compileShaderCodes();
    if (propagate)
      for (int i=childMeshes.length-1; i>-1; i--)
        childMeshes[i].setSpecularMap(img,propagate);
  }//endfunction

  /**
   * brute force check to reduce the number of vertices uploaded by reusing identical vertices
   */
  void compressGeometry([bool propagate=false])
  {
    List<double> oV = _vertData;
    List<int> oI = _idxData;

    if (oV!=null && oI!=null)
    {
      List<VertexData> tV = new List<VertexData>();
      List<int> nI = new List<int>();

      int n = 0;       // tracks length of tV
      int l = oI.length;
      for (int i=0; i<l; i++)
      {
        int oidx = oI[i]*8;    // old index
        double vx = oV[oidx+0];
        double vy = oV[oidx+1];
        double vz = oV[oidx+2];
        double nx = oV[oidx+5];
        double ny = oV[oidx+6];
        double nz = oV[oidx+7];
        double u = oV[oidx+3];
        double v = oV[oidx+4];

        int nidx = -1;      // new index
        n = tV.length;
        for (int j=n-1; j>-1 && nidx==-1; j--)
        {
          VertexData vd = tV[j];
          if (vd.vx==vx && vd.vy==vy && vd.vz==vz &&
              vd.nx==nx && vd.ny==ny && vd.nz==nz &&
              vd.u==u && vd.v==v)
            nidx = j;
        }

        if (nidx==-1)
        {
          nidx = n;
          tV.add(new VertexData(vx,vy,vz,nx,ny,nz,u,v));
        }
        nI.add(nidx);
      }//endfor

      List<double> nV = new List<double>();
      n = tV.length;
      for (int i=0; i<n; i++)
      {
        VertexData vd = tV[i];
        nV.addAll([vd.vx,vd.vy,vd.vz,vd.u,vd.v,vd.nx,vd.ny,vd.nz]);
      }

      setGeometry(nV,nI);
      print("compressed from "+(oV.length~/8).toString()+" to "+(nV.length~/8).toString()+" vertices "+
            (nV.length/oV.length*100).toStringAsFixed(2)+"% of original");
    }

    // ----- do for submeshes too if required
    if (propagate)
    for (int i=0; i<childMeshes.length; i++)
      childMeshes[i].compressGeometry(propagate);
  }//endfunction

  /**
   * creates a texture out of given imageElement
   */
  static Map _uploadedTextures = new Map();
  static Texture createMipMapTexture(ImageElement image,[Texture texture=null])
  {
    if (_uploadedTextures.containsKey(image))
      return _uploadedTextures[image];
    if (_gl==null) return null;
    if (texture==null) texture = _gl.createTexture();
    _gl.pixelStorei(UNPACK_FLIP_Y_WEBGL, 1);
    _gl.bindTexture(TEXTURE_2D, texture);
    _gl.texImage2DImage(TEXTURE_2D, 0, RGBA, RGBA, UNSIGNED_BYTE, image);
    _gl.texParameteri(TEXTURE_2D, TEXTURE_MAG_FILTER, LINEAR);
    _gl.texParameteri(TEXTURE_2D, TEXTURE_MIN_FILTER, LINEAR_MIPMAP_NEAREST);
    _gl.generateMipmap(TEXTURE_2D);
    _uploadedTextures[image] = texture;
    print("#textures = "+_uploadedTextures.length.toString());
    return texture;
  }//endfunction

  /**
   * renders the given mesh tree
   */
  static void render(RenderingContext context,Mesh tree,int viewWidth, int viewHeight, double aspect)
  {
    if (context==null) return;
    _gl = context;
    _gl.clearColor(0.0, 0.0, 0.0, 1.0);  // clear color and alpha

    // Basic viewport setup and clearing of the screen
    _gl.viewport(0, 0, viewWidth, viewHeight);
    _gl.clear(COLOR_BUFFER_BIT | DEPTH_BUFFER_BIT);
    _gl.enable(DEPTH_TEST);
    _gl.enable(BLEND);

    // ----- calculate lighting info
    Float32List lPs = new Float32List(_lightPts.length*3);     // light positions
    Float32List lCs = new Float32List(_lightPts.length*3);     // light colors
    for (int i=_lightPts.length-1; i>-1; i--)
    {
      PointLight lpt = _lightPts[i];
      lPs[i*3+0] = lpt.px*_viewT.aa+lpt.py*_viewT.ab+lpt.pz*_viewT.ac+_viewT.ad;  // transformed lightPt x
      lPs[i*3+1] = lpt.px*_viewT.ba+lpt.py*_viewT.bb+lpt.pz*_viewT.bc+_viewT.bd;  // transformed lightPt y
      lPs[i*3+2] = lpt.px*_viewT.ca+lpt.py*_viewT.cb+lpt.pz*_viewT.cc+_viewT.cd;  // transformed lightPt z
      lCs[i*3+0] = lpt.r;
      lCs[i*3+1] = lpt.g;
      lCs[i*3+2] = lpt.b;
    }//endfor

    List<Mesh> Mshs = new List<Mesh>();
    tree._flattenTree(_viewT,Mshs);

    for (int i=Mshs.length-1; i>-1; i--)
    {
      Mesh m=Mshs[i];   // current mesh to render

      if (m._prog==null && m._vertData!=null && m._idxData!=null)
        m._prepForRender();

      if (m._prog!=null)
      {
        if (m._prepOnState!=_stateId) m._compileShaderCodes();     // recompile shader code if global state changed

        _gl.useProgram(m._prog);                                   // set shader program
        _gl.bindBuffer(ELEMENT_ARRAY_BUFFER, m._idxBuffer);        // bind index buffer
        _gl.bindBuffer(ARRAY_BUFFER, m._vertBuffer);               // bind vertex buffer
        if (m._diffMap!=null)
        {
          _gl.activeTexture(TEXTURE0);
          _gl.bindTexture(TEXTURE_2D, m._diffMap);
          _gl.uniform1i(_gl.getUniformLocation(m._prog, 'diffMap'), 0);
        }
        if (m._normMap!=null)
        {
          _gl.activeTexture(TEXTURE1);
          _gl.bindTexture(TEXTURE_2D, m._normMap);
          _gl.uniform1i(_gl.getUniformLocation(m._prog, 'normMap'), 1);
        }
        if (m._specMap!=null)
        {
          _gl.activeTexture(TEXTURE2);
          _gl.bindTexture(TEXTURE_2D, m._specMap);
          _gl.uniform1i(_gl.getUniformLocation(m._prog, 'specMap'), 2);
        }
        int aLoc = _gl.getAttribLocation(m._prog, 'aPosn');
        _gl.enableVertexAttribArray(aLoc);
        _gl.vertexAttribPointer(aLoc, 3, FLOAT, false, 11*4, 0);    // specify vertex
        aLoc = _gl.getAttribLocation(m._prog, 'aUV');
        _gl.enableVertexAttribArray(aLoc);
        _gl.vertexAttribPointer(aLoc, 2, FLOAT, false, 11*4, 3*4);  // specify UV
        aLoc = _gl.getAttribLocation(m._prog, 'aNorm');
        _gl.enableVertexAttribArray(aLoc);
        _gl.vertexAttribPointer(aLoc, 3, FLOAT, false, 11*4, 5*4);  // specify normal
        aLoc = _gl.getAttribLocation(m._prog, 'aTang');
        _gl.enableVertexAttribArray(aLoc);
        _gl.vertexAttribPointer(aLoc, 3, FLOAT, false, 11*4, 8*4);  // specify tangent

        _gl.uniformMatrix4fv(_gl.getUniformLocation(m._prog, 'uPMatrix'), false, _persT.buf());
        _gl.uniformMatrix4fv(_gl.getUniformLocation(m._prog, 'uMVMatrix'), false, m._workingT.buf());

        _gl.uniform3fv(_gl.getUniformLocation(m._prog, 'lPoints'),lPs);  // upload light posns
        _gl.uniform3fv(_gl.getUniformLocation(m._prog, 'lColors'),lCs);  // upload light colors

        //_gl.drawArrays(TRIANGLES, 0, m.numTris*3);
        _gl.drawElements(TRIANGLES, m.numTris*3, UNSIGNED_SHORT, 0);
      }//endif
    }//endfor

  }//endfunction

  /**
   * compiles and uploads shader program
   */
  static Map _uploadedPrograms = new Map();
  static Program _uploadShaderProgram(String vertSrc, String fragSrc)
  {
    if (_uploadedPrograms.containsKey(vertSrc+"\n"+fragSrc))
      return _uploadedPrograms[vertSrc+"\n"+fragSrc];

    if (_gl==null) return null;

    Shader vertShader = _gl.createShader(VERTEX_SHADER);
    _gl.shaderSource(vertShader, vertSrc);
    _gl.compileShader(vertShader);

    Shader fragShader = _gl.createShader(FRAGMENT_SHADER);
    _gl.shaderSource(fragShader, fragSrc);
    _gl.compileShader(fragShader);

    Program prog = _gl.createProgram();
    _gl.attachShader(prog, vertShader);
    _gl.attachShader(prog, fragShader);
    _gl.linkProgram(prog);

    Object linkStat = _gl.getProgramParameter(prog, LINK_STATUS);
    print("linkStatus=" + linkStat.toString());
    if (linkStat.toString()!='true') print("vertSrc:\n"+vertSrc+"\nfragSrc:\n"+fragSrc);

    _uploadedPrograms[vertSrc+"\n"+fragSrc] = prog;
    print("#programs = " + _uploadedPrograms.length.toString());

    return prog;
  }//endfunction

  /**
   * sets camera to look from (px,py,pz) at (tx,ty,tz), camera is always oriented y up with elevation angle
   */
  static Matrix4 setCamera(double px,double py,double pz,double tx,double ty,double tz,[double fov=45.0,double near=1.0,double far=1000.0])
  {
    double nearClip = max(0,near);
    double farClip = max(0,far);
    nearClip = min(near,far);
    farClip = max(near,far);
    double aspectRatio = 1.0;
    _persT = Matrix4.perspective(fov, aspectRatio, nearClip, farClip);
    _viewT = getViewTransform(px,py,pz,tx,ty,tz);
    _camT = _viewT.inverse();
    return _camT.scale(1.0,1.0,1.0); // duplicate and return
  }//endfunction

  /**
   * calculate tangents for normal mapping, quite heavy calculation
   * input: [vx,vy,vz,u,v,nx,ny,nz,...]
   * output: [vx,vy,vz,u,v,nx,ny,nz,tx,ty,tz,...]
   */
  static List<Vector3> _TBRV=null;
  static List<double> calcTangentBasis(List<int> idxs,List<double> vData)
  {
    /*
    let a be vector from p to q
    let b be vector from p to r

    p(ax,ay) + q(bx,by) s.t    (y axis)
    p*ay + q*by = 1  ... (1)
    p*ax + q*bx = 0  ... (2)

    p*ax = -q*bx
    p = -q*bx/ax   ... (2a)
    sub in (1)

    -q*ay*bx/ax + q*by = 1
    q = 1/(by-ay*bx/ax)
    */

    int i=0;
    int n=vData.length~/8;
    Vector3 v = null;

    if (_TBRV==null) _TBRV=new List<Vector3>();
    for (i=_TBRV.length-1; i>=0; i--)  {v=_TBRV[i]; v.x=0.0; v.y=0.0; v.z=0.0;} // reset vector
    for (i=_TBRV.length; i<n; i++) _TBRV.add(new Vector3(0.0,0.0,0.0));

    n = idxs.length;
    for (i=0; i<n;) // for each triangle
    {
      int i0 = idxs[i];  i++;  // tri point index 0
      int i1 = idxs[i];  i++;  // tri point index 1
      int i2 = idxs[i];  i++;  // tri point index 2

      double pax = vData[i1*8+3] - vData[i0*8+3];
      double ax = pax;
      int tmp;
      do {
        tmp=i0; i0=i1; i1=i2; i2=tmp;
        ax = vData[i1*8+3] - vData[i0*8+3];
      } while (ax*ax>pax*pax);
      tmp=i2; i2=i1; i1=i0; i0=tmp;

      int p0 = i0*8+3;
      int p1 = i1*8+3;
      int p2 = i2*8+3;
      double tx = vData[p0++];
      double ty = vData[p0];
             ax = vData[p1++] - tx; // vector a in uv space
      double ay = vData[p1] - ty;
      double bx = vData[p2++] - tx; // vector b in uv space
      double by = vData[p2] - ty;
      double q = 1/(by-ay*bx/ax);
      double p = -q*bx/ax;

      // find tangent vector from p q
      p0 = i0*8;
      p1 = i1*8;
      p2 = i2*8;
      tx = vData[p0++];
      ty = vData[p0++];
      double tz = vData[p0];
      ax = vData[p1++] - tx;
      ay = vData[p1++] - ty;
      double az = vData[p1] - tz;   // vector a in object space
      p0 = i0*8;
      bx = vData[p2++] - tx;
      by = vData[p2++] - ty;
      double bz = vData[p2] - tz;   // vector b in object space

      tx = p*ax+q*bx;
      ty = p*ay+q*by;
      tz = p*az+q*bz;
      v = _TBRV[i0];   v.x+=tx; v.y+=ty; v.z+=tz; v.w++;
      v = _TBRV[i1];   v.x+=tx; v.y+=ty; v.z+=tz; v.w++;
      v = _TBRV[i2];   v.x+=tx; v.y+=ty; v.z+=tz; v.w++;
    }//endfor

    // ----- get tangent results for each corresponding point
    List<double> R = new List<double>();
    n = vData.length~/8;
    for (i=0; i<n; i++)
    {
      v = _TBRV[i];
      int p0 = i*8+5;
      double ax = vData[p0++];
      double ay = vData[p0++];
      double az = vData[p0];

      double tx = v.y*az - v.z*ay; // cross product tangent
      double ty = v.z*ax - v.x*az;
      double tz = v.x*ay - v.y*ax;
      double tl = 1/sqrt(tx*tx+ty*ty+tz*tz);
      tx*=tl; ty*=tl; tz*=tl;

      p0 = i*8;
      R.addAll([vData[p0++],vData[p0++],vData[p0++],  // vx,vy,vz
                vData[p0++],vData[p0++],              // u,v
                vData[p0++],vData[p0++],vData[p0++],  // nx,ny,nz
                tx,ty,tz]);               // tx,ty,tz
    }//endfor

    return R;
  }//endfunction

  /**
   * sets global point lights positions and colors
   */
  static void setLights(List<PointLight> lights)
  {
    if (_lightPts.length!=lights.length) _stateId++; // to trigger recompile of shader codes
    _lightPts=lights;
  }//endfunction

  /**
   * returns scene transform matrix eqv of camera looking from (px,py,pz) at point (tx,ty,tz)
   */
  static Matrix4 getViewTransform(double px,double py,double pz,double tx,double ty,double tz)
  {
    double vx = tx-px;
    double vy = ty-py;
    double vz = tz-pz;
    double roty = atan2(vx,vz);
    double rotx = atan2(-vy,sqrt(vx*vx+vz*vz));
    return new Matrix4().translate(px,py,pz).rotY(-roty).rotX(-rotx);
  }//endfunction

  /**
   * creates a texture tetra from given square texture
   */
  static Mesh createTetra([double l=1.0,bool soft=true])
  {
    // ----- front point
    double ax = 0.0;
    double ay = l/2/sqrt(3)*2;
    double az = 0.0;

    // ----- back left point
    double bx =-l/2;
    double by =-l/2/sqrt(3);
    double bz = 0.0;

    // ----- back right point
    double cx = l/2;
    double cy =-l/2/sqrt(3);
    double cz = 0.0;

    // ----- top point
    double dx = 0.0;
    double dy = 0.0;
    double dz = sqrt( l*l/4*3 - cy*cy );

    az-=l-dz;
    bz-=l-dz;
    cz-=l-dz;
    dz-=l-dz;

    List<double> VData = new List<double>();
    if (soft)
    {
      VData.addAll([ax,ay,az, 1/2,1-1/2/sqrt(3), ax,ay,az,
                    cx,cy,cz, 0.0,1.0,           cx,cy,cz,
                    bx,by,bz, 1.0,1.0,           bx,by,bz]);
      VData.addAll([dx,dy,dz, 1/2,1-1/2/sqrt(3), dx,dy,dz,
                    ax,ay,az, 0.0,1.0,           ax,ay,az,
                    bx,by,bz, 1.0,1.0,           bx,by,bz]);
      VData.addAll([dx,dy,dz, 1/2,1-1/2/sqrt(3), dx,dy,dz,
                    cx,cy,cz, 0.0,1.0,           cx,cy,cz,
                    ax,ay,az, 1.0,1.0,           ax,ay,az]);
      VData.addAll([dx,dy,dz, 1/2,1-1/2/sqrt(3), dx,dy,dz,
                    bx,by,bz, 0.0,1.0,           bx,by,bz,
                    cx,cy,cz, 1.0,1.0,           cx,cy,cz]);
    }
    else
    {
      VData.addAll([ax,ay,az, 1/2,1-1/2/sqrt(3), 0.0,0.0,0.0,
                    cx,cy,cz, 0.0,1.0,           0.0,0.0,0.0,
                    bx,by,bz, 1.0,1.0,           0.0,0.0,0.0]);
      VData.addAll([dx,dy,dz, 1/2,1-1/2/sqrt(3), 0.0,0.0,0.0,
                    ax,ay,az, 0.0,1.0,           0.0,0.0,0.0,
                    bx,by,bz, 1.0,1.0,           0.0,0.0,0.0]);
      VData.addAll([dx,dy,dz, 1/2,1-1/2/sqrt(3), 0.0,0.0,0.0,
                    cx,cy,cz, 0.0,1.0,           0.0,0.0,0.0,
                    ax,ay,az, 1.0,1.0,           0.0,0.0,0.0]);
      VData.addAll([dx,dy,dz, 1/2,1-1/2/sqrt(3), 0.0,0.0,0.0,
                    bx,by,bz, 0.0,1.0,           0.0,0.0,0.0,
                    cx,cy,cz, 1.0,1.0,           0.0,0.0,0.0]);
    }
    Mesh M = new Mesh(VData);
    if (soft) M.compressGeometry();
    return M;
  }//endfunction

  /**
   * creates a cuboid of specified dimensions
   */
  static Mesh createCube([double w=1.0,double h=1.0,double d=1.0,bool soft=true])
  {
    w/=2.0;
    h/=2.0;
    d/=2.0;
    List<double> V = [-w,-h,-d,  w,-h,-d,  w,h,-d,  -w,h,-d,
                      -w,-h, d,  w,-h, d,  w,h, d,  -w,h, d];

    List<int> I = [0,3,1, 1,3,2,  // front
                   1,2,5, 5,2,6,  // right
                   5,6,4, 4,6,7,  // back
                   4,7,0, 0,7,3,  // left
                   4,0,5, 5,0,1,  // top
                   3,7,2, 2,7,6]; // bottom

    List<double> U = [0.0,1.0, 0.0,0.0, 1.0,1.0, 1.0,1.0, 0.0,0.0, 1.0,0.0];

    int i=0;
    int ul=U.length;
    List<double> VData = new List<double>();
    if (soft)
    {
      for (i=0; i<I.length; i+=3)
      VData.addAll([V[I[i+0]*3+0],V[I[i+0]*3+1],V[I[i+0]*3+2],  // vertex a
                    U[(i*2)%ul+0],U[(i*2)%ul+1],
                    V[I[i+0]*3+0],V[I[i+0]*3+1],V[I[i+0]*3+2],  // normal a
                    V[I[i+1]*3+0],V[I[i+1]*3+1],V[I[i+1]*3+2],  // vertex b
                    U[(i*2)%ul+2],U[(i*2)%ul+3],
                    V[I[i+1]*3+0],V[I[i+1]*3+1],V[I[i+1]*3+2],  // normal b
                    V[I[i+2]*3+0],V[I[i+2]*3+1],V[I[i+2]*3+2],  // vertex c
                    U[(i*2)%ul+4],U[(i*2)%ul+5],
                    V[I[i+2]*3+0],V[I[i+2]*3+1],V[I[i+2]*3+2]]);// normal c
    }
    else
    {
      for (i=0; i<I.length; i+=3)
      VData.addAll([V[I[i+0]*3+0],V[I[i+0]*3+1],V[I[i+0]*3+2],  // vertex a
                    U[(i*2)%ul+0],U[(i*2)%ul+1],
                    0.0,0.0,0.0,  // normal a
                    V[I[i+1]*3+0],V[I[i+1]*3+1],V[I[i+1]*3+2],  // vertex b
                    U[(i*2)%ul+2],U[(i*2)%ul+3],
                    0.0,0.0,0.0,  // normal b
                    V[I[i+2]*3+0],V[I[i+2]*3+1],V[I[i+2]*3+2],  // vertex c
                    U[(i*2)%ul+4],U[(i*2)%ul+5],
                    0.0,0.0,0.0]); // normal c
    }
    Mesh M = new Mesh(VData);
    if (soft) M.compressGeometry();
    return M;
  }//endfunction

  /**
   * creates a sphere of radius r
   */
  static Mesh createSphere(double r,[int lon=32,int lat=16,bool soft=true])
  {
    List<double> S = new List<double>();
    int i=0;
    while (i<lat)
    {
      double sinL0 = sin(PI*i/lat);
      double sinL1 = sin(PI*(i+1)/lat);
      double cosL0 = cos(PI*i/lat);
      double cosL1 = cos(PI*(i+1)/lat);
      List<double> A = createTrianglesBand(sinL0*r,sinL1*r,-cosL0*r,-cosL1*r,lon,soft);

      // ----- adjust UVs of mesh to wrap entire sphere instead
      for (int j=0; j<A.length; j+=8)
      {
        A[j+4]=i/lat+A[j+4]/lat;
        if (soft)
        { // recalculate normals
          double nx = A[j+0];
          double ny = A[j+1];
          double nz = A[j+2];
          double nl = sqrt(nx*nx+ny*ny+nz*nz);
          nx/=nl; ny/=nl; nz/=nl;
          A[j+5]=nx; A[j+6]=ny; A[j+7]=nz;
        }
      }
      S.addAll(A);
      i++;
    }//endfor

    Mesh M = new Mesh(S);
    if (soft) M.compressGeometry();
    return M;
  }//endfunction

  /**
   * create a doughnut shape with band radius r1, thickness r2, of m segments and made of n cylinders
   */
  static Mesh createTorus(double r1,double r2,[int m=32,int n=8,bool soft=true])
  {
    List<double> T = new List<double>();
    int i=0;
    while (i<n)
    {
      List<double> A = createTrianglesBand( r1-r2*cos(i/n*PI*2),
                                            r1-r2*cos((i+1)/n*PI*2),
                                            -r2*sin(i/n*PI*2),
                                            -r2*sin((i+1)/n*PI*2),
                                            m,soft);

      // ----- adjust UVs of mesh to wrap entire torus instead
      List<int> F = [0,0,1,1,0,1];
      for (int j=0; j<A.length; j+=8)
      {
        A[j+7]=i/n+A[j+7]/n;
        if (soft)
        { // 112212
          int stp = F[(j/8).floor()%6]; // step +1 or step
          int lon = (j/48).floor();     // current longitude
          double ang = (lon+stp)/m*PI*2;
          double cx = r1*sin(ang);
          double cy = r1*cos(ang);
          double cz = 0.0;
          double nx = A[j+0]-cx;
          double ny = A[j+1]-cy;
          double nz = A[j+2]-cz;
          double nl = sqrt(nx*nx+ny*ny+nz*nz);
          nx/=nl; ny/=nl; nz/=nl;

          A[j+5]=nx; A[j+6]=ny; A[j+7]=nz;  // override normals
        }
      }
      T.addAll(A);
      i++;
    }//endfunction

    Mesh M = new Mesh(T);
    if (soft) M.compressGeometry();
    return M;
  }//endfunction

  /**
   * creates a circular band of triangles of specified r1,r2 z1,z2
   */
  static List<double> createTrianglesBand(double r1,double r2,double z1,double z2,int n,[bool soft=true])
  {
    List<double> A = new List<double>();
    int i=0;
    while (i<n)
    {
      double a1 = i/n*PI*2;
      double a2 = (i+1)/n*PI*2;

      double sin_a1 = sin(a1);
      double sin_a2 = sin(a2);
      double cos_a1 = cos(a1);
      double cos_a2 = cos(a2);

      if (soft) // apply Smooth Shading
      {
        double mz = (z1+z2)/2;
        if (r2>0) A.addAll([sin_a1*r1,cos_a1*r1,z1, // vertex
                            i/n,0.0,
                            sin_a1*r1,cos_a1*r1,z1, // normal
                            sin_a1*r2,cos_a1*r2,z2, // vertex
                            i/n,1.0,
                            sin_a1*r2,cos_a1*r2,z2, // normal
                            sin_a2*r2,cos_a2*r2,z2, // vertex
                            (i+1)/n,1.0,
                            sin_a2*r2,cos_a2*r2,z2]);// normal
        if (r1>0) A.addAll([sin_a2*r1,cos_a2*r1,z1, // vertex
                            (i+1)/n,0.0,
                            sin_a2*r1,cos_a2*r1,z1, // normal
                            sin_a1*r1,cos_a1*r1,z1, // vertex
                            i/n,0.0,
                            sin_a1*r1,cos_a1*r1,z1, // normal
                            sin_a2*r2,cos_a2*r2,z2, // vertex
                            (i+1)/n,1.0,
                            sin_a2*r2,cos_a2*r2,z2]);// normal
      }
      else
      {
        if (r2>0) A.addAll([sin_a1*r1,cos_a1*r1,z1, // vertex
                            i/n,0.0,
                            0.0,0.0,0.0,            // normal
                            sin_a1*r2,cos_a1*r2,z2, // vertex
                            i/n,1,
                            0.0,0.0,0.0,            // normal
                             sin_a2*r2,cos_a2*r2,z2,// vertex
                            (i+1)/n,1,
                            0.0,0.0,0.0]);          // normal
        if (r1>0) A.addAll([sin_a2*r1,cos_a2*r1,z1, // vertex
                            (i+1)/n,0.0,
                            0.0,0.0,0.0,            // normal
                            sin_a1*r1,cos_a1*r1,z1, // vertex
                            i/n,0.0,
                            0.0,0.0,0.0,            // normal
                            sin_a2*r2,cos_a2*r2,z2, // vertex
                            (i+1)/n,1.0,
                            0.0,0.0,0.0]);          // normal
      }
      i++;
    }//endfor

    return A;
  }//endfunction

  /**
   * Callbacks with parsed obj m as parameter
   */
  static void loadObj(String url,handle(Mesh m))
  {
    List<String> A = url.split(".");
    A.removeLast();
    String mtlUrl = A.join(".")+".mtl";
    A = url.split("/");
    A.removeLast();
    String urlPath = A.join("/")+"/";

    HttpRequest.getString(mtlUrl).then((String mtlData) {
      parseMtlToList(mtlData,urlPath,(List Mtls) {
        HttpRequest.getString(url).then((String objData) {
          handle(Mesh.parseObj(objData,Mtls));
        });
      });
    });
  }//endfunction

  /**
   * parses a given obj format string data s to mesh with null texture
   * Mtls: [id1,bmd1,id2,bmd2,...]
   */
  static Mesh parseObj(String s,[List Mtls=null])
  {
    int i = 0;
    int j = 0;

    // ----- read data from string
    List<String> D = s.split('\n');
    List<double> V = new List<double>();    // list to contain vertices data
    List<double> T = new List<double>();    // list to contain texture coordinates data
    List<double> N = new List<double>();    // list to contain normals data

    List F = new List();              // list to contain triangle faces data
    List G = new List();              // groups array, containing submeshes faces
    List A = null;                    // temp array

    int n = D.length;
    for (i=0; i<n; i++)
    {
      if (D[i].startsWith('v '))       // ----- if position definition
      {
        A = (D[i].substring(2)).split(' ');
        for (j=A.length-1; j>=0; j--)
          if (A[j]=="") A.removeAt(j);
        for (j=0; j<A.length && j<3; j++)
          V.add(double.parse(A[j]));
      }
      else if (D[i].startsWith('vt '))   // ----- if vertex uv definition
      {
        A = (D[i].substring(2)).split(' ');
        for (j=A.length-1; j>=0; j--)
          if (A[j]=="") A.removeAt(j);
        for (j=0; j<A.length && j<2; j++)   // restrict to u,v instead of u,v,t
          T.add(double.parse(A[j]));
      }
      else if (D[i].startsWith('vn '))   // ----- if vertex normal definition
      {
        A = (D[i].substring(2)).split(" ");
        for (j=A.length-1; j>=0; j--)
        {
          //if (A[j].indexOf("e-")!=-1)
            //A[j]=A[j].split("e-")[0];
          if (A[j]=="") A.removeAt(j);
        }
        for (j=0; j<A.length && j<3; j++)
          N.add(double.parse(A[j]));
      }
      else if (D[i].startsWith('f '))    // ----- if face definition
      {
        A = (D[i].substring(2)).split(" ");  // ["v/uv/n","v/uv/n","v/uv/n"]
        for (j=A.length-1; j>=0; j--)
          if (A[j]=="")
            A.removeAt(j);
          else
          {
            while (A[j].split("/").length<3)  A[j] = A[j]+"/-";
            A[j] = A[j].split("//").join("/-/");  // replace null values with "-"
            A[j] = A[j].split("/");   // format of f : [[v,uv,n],[v,uv,n],[v,uv,n]]
          }
        F.add(A);
      }
      else if (D[i].startsWith('o '))    // ----- if object definition
      {
        G.add(F);
        F = new List();
      }
      else if (D[i].startsWith('g '))    // ----- if group definition
      {
        G.add(F);
        F = new List();
      }
      else if (D[i].startsWith('usemtl '))
      {
        F.add((D[i].substring(2)).split(' ')[1]); // material id (defined in mtl file)
      }
    }//endfor

    G.add(F);

    Mesh mmesh = new Mesh();              // main mesh to add all submeshes into

    ImageElement mtl = null;            //
    if (Mtls!=null && Mtls.length>=2) mtl=Mtls[1];  // default material to use to first material

    double onNumParseError(String s) {return 1.0;};

    for (int g=0; g<G.length; g++)
    {
      F = G[g];

      // ----- import faces data -----------------------------
      List<double> verticesData = new List<double>();  // to contain [vx,vy,vz,nx,ny,nz,u,v, ....]
      for (i=0; i<F.length; i++)
      {
        if (F[i] is String) // switch to another material
        {
          if (Mtls!=null) mtl = Mtls[Mtls.indexOf(F[i])+1];
        }
        else
        {
          List f = F[i];        // data of a face: [[v,uv,n],[v,uv,n],[v,uv,n],...]

          for (j=0; j<f.length; j++)
          {
            List p = f[j];      // data of a point: [v,uv,n] to double
            for (int k=0; k<p.length; k++)
              p[k] = (double.parse(p[k],onNumParseError)).toInt()-1;
          }

          // ----- triangulate higher order polygons
          while (f.length>=3)
          {
            A = new List();
            for (j=0; j<3; j++)
              A.addAll(f[j]);
            // A: [v,uv,n,v,uv,n,v,uv,n]

            // ----- get vertices --------------------------------
            double vax = V[A[0]*3+0];
            double vay = V[A[0]*3+1];
            double vaz = V[A[0]*3+2];
            double vbx = V[A[3]*3+0];
            double vby = V[A[3]*3+1];
            double vbz = V[A[3]*3+2];
            double vcx = V[A[6]*3+0];
            double vcy = V[A[6]*3+1];
            double vcz = V[A[6]*3+2];

            // ----- get normals ---------------------------------
            double px = vbx - vax;
            double py = vby - vay;
            double pz = vbz - vaz;

            double qx = vcx - vax;
            double qy = vcy - vay;
            double qz = vcz - vaz;
            // normal by determinant
            double nx = py*qz-pz*qy;  //  unit normal x for the triangle
            double ny = pz*qx-px*qz;  //  unit normal y for the triangle
            double nz = px*qy-py*qx;  //  unit normal z for the triangle

            double nax = nx;   // calculated normals
            double nay = ny;
            double naz = nz;
            double nbx = nx;
            double nby = ny;
            double nbz = nz;
            double ncx = nx;
            double ncy = ny;
            double ncz = nz;
            if (N.length>0)
            {
              nax = N[A[2]*3+0];
              nay = N[A[2]*3+1];
              naz = N[A[2]*3+2];
              nbx = N[A[5]*3+0];
              nby = N[A[5]*3+1];
              nbz = N[A[5]*3+2];
              ncx = N[A[8]*3+0];
              ncy = N[A[8]*3+1];
              ncz = N[A[8]*3+2];
            }

            // ----- get UVs -------------------------------------
            double ua = 0.0;
            double va = 0.0;
            double ub = 1.0;
            double vb = 0.0;
            double uc = 0.0;
            double vc = 1.0;
            if (T.length>0)
            {
              ua = T[A[1]*2+0];
              va = 1-T[A[1]*2+1];
              ub = T[A[4]*2+0];
              vb = 1-T[A[4]*2+1];
              uc = T[A[7]*2+0];
              vc = 1-T[A[7]*2+1];
            }

            verticesData.addAll([ vax,vay,vaz, ua,va, nax,nay,naz,
                                  vbx,vby,vbz, ub,vb, nbx,nby,nbz,
                                  vcx,vcy,vcz, uc,vc, ncx,ncy,ncz]);

            f.removeAt(1);
          }//endwhile
        }//endelse
      }//endfor i

      if (verticesData.length>0)
      {
        Mesh cm = new Mesh(verticesData); //
        cm.compressGeometry(false);
        if (mtl!=null) cm.setTexture(mtl);
        mmesh.addChild(cm);
      }
    }//endfor g

    return mmesh; // returns mesh with submeshes in it
  }//endfunction

  /**
   * returns a list of loaded ImageElements
   */
  static void parseMtlToList(String s,String folder,callBack(List l))
  {
    // ----- create an array of [id1,url1,id2,url2,...] ----------------
    List L = [];
    List<String> D = s.split("\n");
    for (int i=0; i<D.length; i++)
    {
      if (D[i].startsWith("newmtl "))
      {
        if (L.length%2==1) L.add(null);   // if prev id does not have mtl
        L.add(D[i].split(" ").last);      // mtl id
      }
      else if (D[i].startsWith("map_Kd "))
      {
        if (L.length%2==1) L.add(D[i].split(" ").last); // tex url, ensure only one id to one mtl
      }
    }//endfor
    print("parseMtlToList L="+L.toString());

    int idx=0;
    void loadNextImg()
    {
      if (idx>=L.length)
        callBack(L);
      else
      {
        print("parseMtlToList loading "+folder+L[idx+1]);
        ImageElement img = new ImageElement(src:folder+L[idx+1]);
        img.onLoad.listen((dat) {
          L[idx+1]=img;
          idx+=2;
          loadNextImg();
        });
      }
    }//endfunction
    loadNextImg();
  }//endfunction
}//endclass

/**
 * data class for poing light info
 */
class PointLight
{
  double px;
  double py;
  double pz;
  double att;
  double r;
  double g;
  double b;

  PointLight(double x,double y,double z,[double red=1.0,double green=1.0,double blue=1.0,double attenuation=100.0])
  {
    px=x;
    py=y;
    pz=z;
    r=red;
    g=green;
    b=blue;
    att=attenuation;
  }//endconstr
}//endclass

/**
 * convenience class to store vertices info
 */
class VertexData
{
  double vx; // for vertex
  double vy;
  double vz;
  double nx; // for normal
  double ny;
  double nz;
  double tx; // for tangent
  double ty;
  double tz;
  double u;  // for UV
  double v;
  double w;
  int idx;   // for weight idx

  VertexData([double vx=0.0,double vy=0.0,double vz=0.0,double nx=0.0,double ny=0.0,double nz=0.0,double u=0.0,double v=0.0,double w=0.0,int idx=0,double tx=0.0,double ty=0.0,double tz=0.0])
  {
    this.vx = vx;
    this.vy = vy;
    this.vz = vz;
    this.nx = nx;
    this.ny = ny;
    this.nz = nz;
    this.tx = tx;
    this.ty = ty;
    this.tz = tz;
    this.u = u;
    this.v = v;
    this.w = w;
    this.idx = idx;
  }//endconstructor

  VertexData clone()
  {
    return new VertexData(vx,vy,vz,nx,ny,nz,u,v,w,idx,tx,ty,tz);
  }//endfunction

  String toString()
  {
    return  vx.toStringAsFixed(2)+","+vy.toStringAsFixed(2)+","+vz.toStringAsFixed(2)+","+
            nx.toStringAsFixed(2)+","+ny.toStringAsFixed(2)+","+nz.toStringAsFixed(2)+","+
            u.toStringAsFixed(2)+","+v.toStringAsFixed(2)+","+w.toStringAsFixed(2)+","+idx.toString();
  }
}//endclass


// Trimmed down matrix code - its best you use a well known library
// instead of this code.

/**
 * 3 dimensional vector.
 */
class Vector3
{
  double x,y,z,w;

  Vector3([double x_=0.0, double y_=0.0, double z_=0.0])
  {
    x = x_;
    y = y_;
    z = z_;
    w = 0.0;
  }

  double length() => sqrt(x*x + y*y + z*z);

  Vector3 normalize()
  {
    double len = length();
    if (len!=0.0)
      return new Vector3(x/len, y/len, z/len);
    return new Vector3(x, y, z);
  }

  Vector3 operator -()
  {
    return new Vector3(-x, -y, -z);
  }

  Vector3 operator -(Vector3 other)
  {
    return new Vector3(x - other.x, y - other.y, z - other.z);
  }

  Vector3 cross(Vector3 other)
  {
    double nx = y * other.z - z * other.y;
    double ny = z * other.x - x * other.z;
    double nz = x * other.y - y * other.x;
    return new Vector3(nx, ny, nz);
  }

  Vector3 scale(double by)
  {
    return new Vector3(x*by, y*by, z*by);
  }

  String toString()
  {
    return "($x,$y,$z)";
  }
}

/**
* 4x4 matrix class, Transform matrix for Mesh class
*/
class Matrix4
{
  double aa,ab,ac,ad;
  double ba,bb,bc,bd;
  double ca,cb,cc,cd;
  double da,db,dc,dd;

  Matrix4([double aa_=1.0,double ab_=0.0,double ac_=0.0,double ad_=0.0,
            double ba_=0.0,double bb_=1.0,double bc_=0.0,double bd_=0.0,
            double ca_=0.0,double cb_=0.0,double cc_=1.0,double cd_=0.0,
            double da_=0.0,double db_=0.0,double dc_=0.0,double dd_=1.0])
  {
    aa = aa_;  ab = ab_;  ac = ac_;  ad = ad_;
    ba = ba_;  bb = bb_;  bc = bc_;  bd = bd_;
    ca = ca_;  cb = cb_;  cc = cc_;  cd = cd_;
    da = da_;  db = db_;  dc = dc_;  dd = dd_;
  }//end constructor

  /**
  * returns if this matrix equals given matrix
  */
  bool equals(Matrix4 M)
  {
    return  aa==M.aa && ab==M.ab && ac==M.ac && ad==M.ad &&
            ba==M.ba && bb==M.bb && bc==M.bc && bd==M.bd &&
            ca==M.ca && cb==M.cb && cc==M.cc && cd==M.cd &&
            da==M.da && db==M.db && dc==M.dc && dd==M.dd;
  }//endfunction

  /**
  * this=A, M=B, mult -> AB
  */
  Matrix4 mult(Matrix4 M)
  {
    double naa = aa*M.aa + ab*M.ba + ac*M.ca + ad*M.da;
    double nab = aa*M.ab + ab*M.bb + ac*M.cb + ad*M.db;
    double nac = aa*M.ac + ab*M.bc + ac*M.cc + ad*M.dc;
    double nad = aa*M.ad + ab*M.bd + ac*M.cd + ad*M.dd;

    double nba = ba*M.aa + bb*M.ba + bc*M.ca + bd*M.da;
    double nbb = ba*M.ab + bb*M.bb + bc*M.cb + bd*M.db;
    double nbc = ba*M.ac + bb*M.bc + bc*M.cc + bd*M.dc;
    double nbd = ba*M.ad + bb*M.bd + bc*M.cd + bd*M.dd;

    double nca = ca*M.aa + cb*M.ba + cc*M.ca + cd*M.da;
    double ncb = ca*M.ab + cb*M.bb + cc*M.cb + cd*M.db;
    double ncc = ca*M.ac + cb*M.bc + cc*M.cc + cd*M.dc;
    double ncd = ca*M.ad + cb*M.bd + cc*M.cd + cd*M.dd;

    double nda = da*M.aa + db*M.ba + dc*M.ca + dd*M.da;
    double ndb = da*M.ab + db*M.bb + dc*M.cb + dd*M.db;
    double ndc = da*M.ac + db*M.bc + dc*M.cc + dd*M.dc;
    double ndd = da*M.ad + db*M.bd + dc*M.cd + dd*M.dd;

    return new Matrix4(naa,nab,nac,nad,
                        nba,nbb,nbc,nbd,
                        nca,ncb,ncc,ncd,
                        nda,ndb,ndc,ndd);
  }//end Function

  /**
  * append multiplication to target matrix M -> AM
  */
  static Matrix4 appendMult(Matrix4 M,
                              double aa,double ab,double ac,double ad,
                              double ba,double bb,double bc,double bd,
                              double ca,double cb,double cc,double cd,
                              double da,double db,double dc,double dd)
  {
    double naa = aa*M.aa + ab*M.ba + ac*M.ca + ad*M.da;
    double nab = aa*M.ab + ab*M.bb + ac*M.cb + ad*M.db;
    double nac = aa*M.ac + ab*M.bc + ac*M.cc + ad*M.dc;
    double nad = aa*M.ad + ab*M.bd + ac*M.cd + ad*M.dd;

    double nba = ba*M.aa + bb*M.ba + bc*M.ca + bd*M.da;
    double nbb = ba*M.ab + bb*M.bb + bc*M.cb + bd*M.db;
    double nbc = ba*M.ac + bb*M.bc + bc*M.cc + bd*M.dc;
    double nbd = ba*M.ad + bb*M.bd + bc*M.cd + bd*M.dd;

    double nca = ca*M.aa + cb*M.ba + cc*M.ca + cd*M.da;
    double ncb = ca*M.ab + cb*M.bb + cc*M.cb + cd*M.db;
    double ncc = ca*M.ac + cb*M.bc + cc*M.cc + cd*M.dc;
    double ncd = ca*M.ad + cb*M.bd + cc*M.cd + cd*M.dd;

    double nda = da*M.aa + db*M.ba + dc*M.ca + dd*M.da;
    double ndb = da*M.ab + db*M.bb + dc*M.cb + dd*M.db;
    double ndc = da*M.ac + db*M.bc + dc*M.cc + dd*M.dc;
    double ndd = da*M.ad + db*M.bd + dc*M.cd + dd*M.dd;

    M.aa=naa; M.ab=nab; M.ac=nac; M.ad=nad;
    M.ba=nba; M.bb=nbb; M.bc=nbc; M.bd=nbd;
    M.ca=nca; M.cb=ncb; M.cc=ncc; M.cd=ncd;
    M.da=nda; M.db=ndb; M.dc=ndc; M.dd=ndd;

    return M;
  }//endfunction

  /**
  * returns the determinant of this 4x4 matrix
  */
  double determinant()
  {
    return aa*bb*cc*dd + aa*bc*cd*db + aa*bd*cb*db
        + ab*ba*cd*dc + ab*bc*ca*dd + ab*bd*cc*da
        + ac*ba*cb*dd + ac*bb*cd*da + ac*bd*ca*db
        + ad*ba*cc*db + ad*bb*ca*dc + ad*bc*cb*da
        - aa*bb*cd*dc - aa*bc*cb*dd - aa*bd*cc*db
        - ab*ba*cc*dd - ab*bc*cd*da - ab*bd*ca*dc
        - ac*ba*cd*db - ac*bb*ca*dd - ac*bd*cb*da
        - ad*ba*cb*dc - ad*bb*cc*da - ad*bc*ca*db;
  }//endfunction

  /**
  * returns the determinant if the inner 3x3 matrix, which is also the scaling factor
  */
  double determinant3()
  {
    // aei+bfg+cdh-ceg-bdi-afh
    return aa*bb*cc + ab*bc*ca + ac*ba*cb - ac*bb*ca - ab*ba*cc - aa*bc*cb;
  }//endfunction

  /**
  * returns the inverse matrix of this matrix
  */
  Matrix4 inverse()
  {
    double _det = determinant();
    if (_det==0)  return null;
    double naa = bb*cc*dd + bc*cd*db + bd*cb*dc - bb*cd*dc - bc*cb*dd - bd*cc*db;
    double nab = ab*cd*dc + ac*cb*dd + ad*cc*db - ab*cc*dd - ac*cd*db - ad*cb*dc;
    double nac = ab*bc*dd + ac*bd*db + ad*bb*dc - ab*bd*dc - ac*bb*dd - ad*bc*db;
    double nad = ab*bd*cc + ac*bb*cd + ad*bc*cb - ab*bc*cd - ac*bd*cb - ad*bb*cc;
    double nba = ba*cd*dc + bc*ca*dd + bd*cc*da - ba*cc*dd - bc*cd*da - bd*ca*dc;
    double nbb = aa*cc*dd + ac*cd*da + ad*ca*dc - aa*cd*dc - ac*ca*dd - ad*cc*da;
    double nbc = aa*bd*dc + ac*ba*dd + ad*bc*da - aa*bc*dd - ac*bd*da - ad*ba*dc;
    double nbd = aa*bc*cd + ac*bd*ca + ad*ba*cc - aa*bd*cc - ac*ba*cd - ad*bc*ca;
    double nca = ba*cb*dd + bb*cd*da + bd*ca*db - ba*cd*db - bb*ca*dd - bd*cb*da;
    double ncb = aa*cd*db + ab*ca*dd + ad*cb*da - aa*cb*dd - ab*cd*da - ad*ca*db;
    double ncc = aa*bb*dd + ab*bd*da + ad*ba*db - aa*bd*db - ab*ba*dd - ad*bb*da;
    double ncd = aa*bd*cb + ab*ba*cd + ad*bb*ca - aa*bb*cd - ab*bd*ca - ad*ba*cb;
    double nda = ba*cc*db + bb*ca*dc + bc*cb*da - ba*cb*dc - bb*cc*da - bc*ca*db;
    double ndb = aa*cb*dc + ab*cc*da + ac*ca*db - aa*cc*db - ab*ca*dc - ac*cb*da;
    double ndc = aa*bc*db + ab*ba*dc + ac*bb*da - aa*bb*dc - ab*bc*da - ac*ba*db;
    double ndd = aa*bb*cc + ab*bc*ca + ac*ba*cb - aa*bc*cb - ab*ba*cc - ac*bb*ca;
    _det = 1/_det;  // determinant inverse, to prevent 16 divisions
    return new Matrix4( naa*_det,nab*_det,nac*_det,nad*_det,
                        nba*_det,nbb*_det,nbc*_det,nbd*_det,
                        nca*_det,ncb*_det,ncc*_det,ncd*_det,
                        nda*_det,ndb*_det,ndc*_det,ndd*_det);
  }//endfunction

  /*
  * returns new matrix with exact cloned values
  */
  Matrix4 clone()
  {
    return new Matrix4(aa,ab,ac,ad, ba,bb,bc,bd, ca,cb,cc,cd, da,db,dc,dd);
  }//endfunction

  /**
  * returns new transform matrix scaled by (xs,ys,zs)
  */
  Matrix4 scale([double xs=1.0,double ys=1.0,double zs=1.0])
  {
    double naa = xs;
    double nab = 0.0;
    double nac = 0.0;
    double nad = 0.0;

    double nba = 0.0;
    double nbb = ys;
    double nbc = 0.0;
    double nbd = 0.0;

    double nca = 0.0;
    double ncb = 0.0;
    double ncc = zs;
    double ncd = 0.0;

    double nda = 0.0;
    double ndb = 0.0;
    double ndc = 0.0;
    double ndd = 1.0;

    return appendMult(this.clone(),
                      naa,nab,nac,nad,
                      nba,nbb,nbc,nbd,
                      nca,ncb,ncc,ncd,
                      nda,ndb,ndc,ndd);
  }//end Function

  /**
  * returns new transform matrix rotated about Z
  */
  Matrix4 rotZ(double a)
  {
    double cosA = cos(a);
    double sinA = sin(a);

    double naa = cosA;
    double nab =-sinA;
    double nac = 0.0;
    double nad = 0.0;

    double nba = sinA;
    double nbb = cosA;
    double nbc = 0.0;
    double nbd = 0.0;

    double nca = 0.0;
    double ncb = 0.0;
    double ncc = 1.0;
    double ncd = 0.0;

    double nda = 0.0;
    double ndb = 0.0;
    double ndc = 0.0;
    double ndd = 1.0;

    return appendMult(  this.clone(),
              naa,nab,nac,nad,
              nba,nbb,nbc,nbd,
              nca,ncb,ncc,ncd,
              nda,ndb,ndc,ndd);
  }//end Function rotZ

  /**
  * returns new transform matrix rotated about Y
  */
  Matrix4 rotY(double a)
  {
    double cosA = cos(a);
    double sinA = sin(a);

    double naa = cosA;
    double nab = 0.0;
    double nac = sinA;
    double nad = 0.0;

    double nba = 0.0;
    double nbb = 1.0;
    double nbc = 0.0;
    double nbd = 0.0;

    double nca =-sinA;
    double ncb = 0.0;
    double ncc = cosA;
    double ncd = 0.0;

    double nda = 0.0;
    double ndb = 0.0;
    double ndc = 0.0;
    double ndd = 1.0;

    return appendMult(  this.clone(),
              naa,nab,nac,nad,
              nba,nbb,nbc,nbd,
              nca,ncb,ncc,ncd,
              nda,ndb,ndc,ndd);
  }//end Function rotY

  /**
  * returns new transform matrix rotated about X
  */
  Matrix4 rotX(double a)
  {
    double cosA = cos(a);
    double sinA = sin(a);

    double naa = 1.0;
    double nab = 0.0;
    double nac = 0.0;
    double nad = 0.0;

    double nba = 0.0;
    double nbb = cosA;
    double nbc =-sinA;
    double nbd = 0.0;

    double nca = 0.0;
    double ncb = sinA;
    double ncc = cosA;
    double ncd = 0.0;

    double nda = 0.0;
    double ndb = 0.0;
    double ndc = 0.0;
    double ndd = 1.0;

    return appendMult(this.clone(),
              naa,nab,nac,nad,
              nba,nbb,nbc,nbd,
              nca,ncb,ncc,ncd,
              nda,ndb,ndc,ndd);
  }//end Function rotX

  /**
  * returns new transform matrix rotated by angle specified by 2 vectors
  */
  Matrix4 rotFromTo(double ax,double ay,double az,double bx,double by,double bz)
  {
    double _al = 1/sqrt(ax*ax+ay*ay+az*az);
    ax*=_al;
    ay*=_al;
    az*=_al;
    double _bl = 1/sqrt(bx*bx+by*by+bz*bz);
    bx*=_bl;
    by*=_bl;
    bz*=_bl;

    // ----- reversed direction special case
    if ((ax+bx)*(ax+bx) + (ay+by)*(ay+by) + (az+bz)*(az+bz)<0.0000001)
    {
      Matrix4 fM = new Matrix4(-1.0, 0.0, 0.0, 0.0,
                                  0.0,-1.0, 0.0, 0.0,
                                  0.0, 0.0, 1.0, 0.0,
                                  0.0, 0.0, 0.0, 1.0); // flip from up to down
      if (ay>0)
        fM = fM.mult(new Matrix4().rotFromTo(ax,ay,az,0.0,1.0,0.0)).rotFromTo(0.0,1.0,0.0,ax,ay,az);
      else
        fM = fM.mult(new Matrix4().rotFromTo(ax,ay,az,0.0,-1.0,0.0)).rotFromTo(0.0,-1.0,0.0,ax,ay,az);

      return appendMult(this.clone(),
                fM.aa,fM.ab,fM.ac,fM.ad,
                fM.ba,fM.bb,fM.bc,fM.bd,
                fM.ca,fM.cb,fM.cc,fM.cd,
                fM.da,fM.db,fM.dc,fM.dd);
    }
    // ----- no rotation special case
    else if ((ax-bx)*(ax-bx) + (ay-by)*(ay-by) + (az-bz)*(az-bz)<0.0000001)
    {
      return this.clone();
    }

    // normal by determinant Tn
    double nx = ay*bz-az*by;  //  normal x for the triangle
    double ny = az*bx-ax*bz;  //  normal y for the triangle
    double nz = ax*by-ay*bx;  //  normal z for the triangle
    return rotAbout(nx,ny,nz,acos(ax*bx+ay*by+az*bz));
  }//endfunction

  /**
  * returns new transform matrix rotated about given axis (ux,uy,uz)
  */
  Matrix4 rotAbout(double ux,double uy,double uz,double a)
  {
    double ul = sqrt(ux*ux + uy*uy + uz*uz);
    if (ul==0)  return this;
    ux/=ul;
    uy/=ul;
    uz/=ul;

    double cosA = cos(a);
    double sinA = sin(a);

    double naa = ux*ux + (1-ux*ux)*cosA;
    double nab = ux*uy*(1-cosA) - uz*sinA;
    double nac = ux*uz*(1-cosA) + uy*sinA;
    double nad = 0.0;

    double nba = ux*uy*(1-cosA) + uz*sinA;
    double nbb = uy*uy + (1-uy*uy)*cosA;
    double nbc = uy*uz*(1-cosA) - ux*sinA;
    double nbd = 0.0;

    double nca = ux*uz*(1-cosA) - uy*sinA;
    double ncb = uy*uz*(1-cosA) + ux*sinA;
    double ncc = uz*uz + (1-uz*uz)*cosA;
    double ncd = 0.0;

    double nda = 0.0;
    double ndb = 0.0;
    double ndc = 0.0;
    double ndd = 1.0;

    return appendMult(this.clone(),
                      naa,nab,nac,nad,
                      nba,nbb,nbc,nbd,
                      nca,ncb,ncc,ncd,
                      nda,ndb,ndc,ndd);
  }//end Function rotAbout

  /**
  * given rotation vector, apply rotation to matrix,
  */
  Matrix4 rotate(double rx,double ry,double rz)
  {
    double rl = sqrt(rx*rx+ry*ry+rz*rz);
    return rotAbout(rx,ry,rz,rl);
  }//endfunction

  /**
  * returns new transform matrix translated (tx,ty,tz)
  */
  Matrix4 translate(double tx,double ty,double tz)
  {
    Matrix4 M = this.clone();
    M.ad += tx;
    M.bd += ty;
    M.cd += tz;
    return M;
  }//end Function translate

  /**
  * returns new transformed Vector3
  */
  Vector3 transform(Vector3 v)
  {
    return new Vector3(v.x*aa+v.y*ab+v.z*ac+ad,
                        v.x*ba+v.y*bb+v.z*bc+bd,
                        v.x*ca+v.y*cb+v.z*cc+cd);
  }//endfunction

  /*
  * returns new rotated Vector3
  */
  Vector3 rotateVector(Vector3 v)
  {
    return new Vector3(v.x*aa+v.y*ab+v.z*ac,
                        v.x*ba+v.y*bb+v.z*bc,
                        v.x*ca+v.y*cb+v.z*cc);
  }//endfunction

  /**
  * returns the string printout of the values
  */
  String toString()
  {
    return  "\n"+
        "|"+padN(aa,5)+","+padN(ab,5)+","+padN(ac,5)+","+padN(ad,5)+"|\n"+
        "|"+padN(ba,5)+","+padN(bb,5)+","+padN(bc,5)+","+padN(bd,5)+"|\n"+
        "|"+padN(ca,5)+","+padN(cb,5)+","+padN(cc,5)+","+padN(cd,5)+"|\n"+
        "|"+padN(da,5)+","+padN(db,5)+","+padN(dc,5)+","+padN(dd,5)+"|\n";
  }//endfunction

  /**
  * returns the quaternion representation of the matrix rotation component
  * code from http://www.cs.princeton.edu/~gewang/projects/darth/stuff/quat_faq.html#Q55
  */
  Vector3 rotationQuaternion()
  {
    double t = 1+aa+bb+cc;  // trace of the matrix
    double s = 0.0;
    Vector3 quat = new Vector3();
    if (t>0.000001)
    {
      s = sqrt(t)*2;
      quat.x = (cb-bc)/s;
      quat.y = (ac-ca)/s;
      quat.z = (ba-ab)/s;
      quat.w = 0.25*s;
    }
    else if (t==0)    // added to attempt a fix on 180 rotations...
    {

      if (aa==1)      quat.x=1.0; // rot 180 about x axis
      else if (bb==1) quat.y=1.0; // rot 180 about y axis
      else            quat.z=1.0; // rot 180 about z axis
    }
    else if (aa>bb && aa>cc)  // column 0
    {
      t = 1+aa-bb-cc;
      s = sqrt(t)*2;
      quat.x = 0.25*s;
      quat.y = (ba+ab)/s;
      quat.z = (ac+ca)/s;
      quat.w = (cb-bc)/s;
    }
    else if (bb>cc)
    {
      t = 1+bb-aa-cc;
      s = sqrt(t)*2;
      quat.x = (ba+ab)/s;
      quat.y = 0.25*s;
      quat.z = (cb+bc)/s;
      quat.w = (ac-ca)/s;
    }
    else
    {
      t = 1+cc-aa-bb;
      s = sqrt(t)*2;
      quat.x = (ac+ca)/s;
      quat.y = (cb+bc)/s;
      quat.z = 0.25*s;
      quat.w = (ba-ab)/s;
    }

    return quat;
  }//endfunction

  /**
   * return in list form for buffer upload
   */
  Float32List buf()
  { // weirdly openGl expects marices to be in colmn major
    Float32List l = new Float32List(16);
    l[0]=aa;  l[4]=ab;  l[8]=ac;  l[12]=ad;
    l[1]=ba;  l[5]=bb;  l[9]=bc;  l[13]=bd;
    l[2]=ca;  l[6]=cb;  l[10]=cc; l[14]=cd;
    l[3]=da;  l[7]=db;  l[11]=dc; l[15]=dd;
    return l;
  }

  /**
  * returns the matrix representation of quaternion rotation w + xi + yj + zk
  */
  static Matrix4 quaternionToMatrix(w,x,y,z)
  {
    double l = sqrt(w*w + x*x + y*y + z*z);
    if (l==0.0) return new Matrix4();
    w/=l; x/=l; y/=l; z/=l;
    double naa = 1 - 2*y*y - 2*z*z;
    double nab = 2*x*y - 2*w*z;
    double nac = 2*x*z + 2*w*y;

    double nba = 2*x*y + 2*w*z;
    double nbb = 1 - 2*x*x - 2*z*z;
    double nbc = 2*y*z - 2*w*x;

    double nca = 2*x*z - 2*w*y;
    double ncb = 2*y*z + 2*w*x;
    double ncc = 1 - 2*x*x - 2*y*y;

    return new Matrix4(naa,nab,nac,0.0,
                        nba,nbb,nbc,0.0,
                        nca,ncb,ncc,0.0);
  }//endfunction

  String padN(double v,int n)
  {
    v = (v*100).round()/100;
    String s = v.toString()+"";
    while (s.length<n)  s = " "+s;
    return s;
  }

  /**
   * Makse a 4x4 matrix perspective projection matrix given a field of view and
   * aspect ratio.
   *
   * [fovyDegrees] field of view (in degrees) of the y-axis
   * [aspectRatio] width to height aspect ratio.
   * [zNear] distance to the near clipping plane.
   * [zFar] distance to the far clipping plane.
   */
  static Matrix4 perspective(double fovyDegrees, double aspectRatio, double zNear, double zFar)
  {
    double height = tan(fovyDegrees/180*PI * 0.5) * zNear.toDouble();
    double width = height * aspectRatio.toDouble();

    Matrix4 dest = new Matrix4();
    double two_near = 2.0 * zNear;
    double right_minus_left = width*2;
    double top_minus_bottom = height*2;
    double far_minus_near = zFar - zNear;
    dest.aa = two_near / right_minus_left;
    dest.bb = two_near / top_minus_bottom;
    dest.cc = -(zFar + zNear) / far_minus_near;
    dest.dc = -1.0;
    dest.cd = -(two_near * zFar) / far_minus_near;
    return dest;
  }//endfunction

}//end class
