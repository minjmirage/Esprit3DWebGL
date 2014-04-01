import 'dart:html';
import 'dart:web_gl';
import 'dart:typed_data';
import 'dart:math';
import 'matrix4.dart';

CanvasElement canvas = querySelector("#glCanvas");
RenderingContext gl;
Matrix4 pMatrix;
Mesh scene;

void main() 
{
  // Nab the context we'll be drawing to.
  print("canvas:" + canvas.width.toString() + "x" + canvas.height.toString());
  gl = canvas.getContext3d();
  if (gl == null) return;

  scene = new Mesh();
  
  Mesh itm = Mesh.createTetra();
  scene.addChild(itm);
  itm = Mesh.createTorus(0.5,0.2);
  itm.transform = itm.transform.translate(sin(0), 0.0, cos(0));
  scene.addChild(itm);
  itm = Mesh.createSphere(0.5);
  itm.transform = itm.transform.translate(sin(PI*2/3), 0.0, cos(PI*2/3));
  scene.addChild(itm);
  itm = Mesh.createCube(1.0,1.0,1.0);
  itm.transform = itm.transform.translate(sin(PI*4/3), 0.0, cos(PI*4/3));
  scene.addChild(itm);
  tick(0.0);
}//endmain

tick(double time) 
{
  window.animationFrame.then(tick);
  print("time:"+time.toString());
  double r = 10.0+5*sin(time/500);
  double sinT = sin(sin(time/1000));
  
  Mesh.setLights([new PointLight(10.0*sin(time/500),0.0,10.0*cos(time/500),1.0,0.0,0.0),
                  new PointLight(10.0*sin(time/500+PI*2/3),0.0,10.0*cos(time/500+PI*2/3),0.0,1.0,0.0),
                  new PointLight(10.0*sin(time/500+PI*4/3),0.0,10.0*cos(time/500+PI*4/3),0.0,0.0,1.0)]);
  Mesh.setCamera(r*cos(sinT)*sin(time/1000), r*sin(sinT), r*cos(sinT)*cos(time/1000), 0.0,0.0,0.0);
  
  Mesh.render(gl,scene,800, 600, 800 / 600);
}


/**
 * 
 */
class Mesh 
{
  Program prog;           // the render program for this mesh
  Buffer vertBuffer;      // uploaded vertices data to GPU
  Buffer idxBuffer;       // uploaded indices data to GPU
  List<Mesh> childMeshes; // list of children
  Matrix4 transform;      // transform for this mesh
  Matrix4 _workingT;      // used for rendering 
  int numTris;            // num of triangles to render
  int prepOnState=-1;     // current state id shader programs are compiled on
  
  static int stateId=0;             // incremented when user changes render parameters
  static List<PointLight> lightPts; // 
  static Matrix4 viewT;             // view transform of the whole scene to be rendered
  static Matrix4 camT;              // camera transform, inverse of viewT
  static Matrix4 persT;             // perspective transform matrix
  
  /**
   * creates new mesh from given geometry data
   */
  Mesh([List<double> vertData=null,List<int> idxData=null])
  {
    // ----- set default global states if not initialized
    if (viewT==null || persT==null) Mesh.setCamera(0.0,0.0,10.0,0.0,0.0,0.0);
    if (lightPts==null) lightPts = [new PointLight(0.0,0.0,10.0)];
    
    // ----- set default for this mesh
    transform = new Matrix4();
    childMeshes = new List<Mesh>();
    
    // ----- upload vertex data to buffer
    if (vertData!=null)
    {
      vertBuffer = gl.createBuffer();
      gl.bindBuffer(ARRAY_BUFFER, vertBuffer);
      gl.bufferDataTyped(ARRAY_BUFFER,new Float32List.fromList(vertData),STATIC_DRAW);
      numTris = vertData.length~/24;
      
      if (idxData==null) idxData = new List<int>();
      for (int i=0; i<numTris*3; i++)  idxData.add(i);
      idxBuffer = gl.createBuffer();
      gl.bindBuffer(ELEMENT_ARRAY_BUFFER, idxBuffer);
      gl.bufferDataTyped(ELEMENT_ARRAY_BUFFER,new Uint16List.fromList(idxData),STATIC_DRAW);
    }//endif
    
    prepareForRender();
  }//endconstr
  
  /**
   * compiles the shader program specific to the mesh and upload
   */
  void prepareForRender()
  {
    // ----- construct vertex shader source -------------
    String vertSrc = '''
     attribute vec3 aPosn;
     attribute vec2 aUV;
     attribute vec3 aNorm;
     
     varying vec3 vPosn;
     varying vec2 vUV;
     varying vec3 vNorm;
     
     uniform mat4 uMVMatrix;
     uniform mat4 uPMatrix;
 
     void main(void) {
         vPosn = (uMVMatrix * vec4(aPosn, 1.0)).xyz;
         gl_Position = uPMatrix * vec4(vPosn, 1.0);
         vUV = aUV;
         vNorm = (uMVMatrix * vec4(aNorm, 0.0)).xyz;
     }
    ''';
    
    // ----- construct fragment shader source -----------
    String fragSrc = '''
     precision highp float;
     
     varying vec2 vUV;
     varying vec3 vNorm;
     varying vec3 vPosn;
     
     uniform vec3 lPoints['''+lightPts.length.toString()+'''];
     uniform vec3 lColors['''+lightPts.length.toString()+'''];
     
     void main(void) 
     {
         vec3 norm = normalize(vNorm);
         vec3 diffuse = vec3(0.0,0.0,0.0);
         vec3 texColor;
         vec3 lightDir;
      ''';
    
    for (int i=0; i<lightPts.length; i++)
    fragSrc += '''
         texColor = lColors['''+i.toString()+'''];
         lightDir = normalize(lPoints['''+i.toString()+'''] - vPosn);
         diffuse = diffuse+texColor*dot(lightDir,norm);
      ''';
    
    fragSrc += '''
         gl_FragColor = vec4(diffuse,1.0);
    }''';
    //print(fragSrc);
    prog = uploadShaderProgram(vertSrc,fragSrc);
    
    prepOnState = stateId;  // specify shader code is built based on this state
    
    /*  
          '''
          precision mediump float;
          
          varying vec2 vUV;
          varying vec3 vNorm;
          varying vec3 vPosn;
          
          uniform float uMaterialShininess;
          
          uniform bool uShowSpecularHighlights;
          uniform bool uUseLighting;
          uniform bool uUseTextures;
          
          uniform vec3 uAmbientColor;
          
          uniform vec3 uPointLightingLocation;
          uniform vec3 uPointLightingSpecularColor;
          uniform vec3 uPointLightingDiffuseColor;
          
          uniform sampler2D uSampler;
          
          
          void main(void) {
              vec3 lightWeighting;
              if (!uUseLighting) {
                  lightWeighting = vec3(1.0, 1.0, 1.0);
              } else {
                  vec3 lightDirection = normalize(uPointLightingLocation - vPosn.xyz);
                  vec3 normal = normalize(vNorm);
          
                  float specularLightWeighting = 0.0;
                  if (uShowSpecularHighlights) {
                      vec3 eyeDirection = normalize(-vPosn.xyz);
                      vec3 reflectionDirection = reflect(-lightDirection, normal);
          
                      specularLightWeighting = pow(max(dot(reflectionDirection, eyeDirection), 0.0), uMaterialShininess);
                  }
          
                  float diffuseLightWeighting = max(dot(normal, lightDirection), 0.0);
                  lightWeighting = uAmbientColor
                      + uPointLightingSpecularColor * specularLightWeighting
                      + uPointLightingDiffuseColor * diffuseLightWeighting;
              }
          
              vec4 fragmentColor;
              if (uUseTextures) {
                  fragmentColor = texture2D(uSampler, vec2(vUV.s, vUV.t));
              } else {
                  fragmentColor = vec4(1.0, 1.0, 1.0, 1.0);
              }
              gl_FragColor = vec4(fragmentColor.rgb * lightWeighting, fragmentColor.a);
          }
        ''';
        */
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
  void flattenTree(Matrix4 T,List<Mesh> L)
  {
    _workingT = T.mult(transform); // working transform of this mesh
    L.add(this);
    for (int i=childMeshes.length-1; i>-1; i--)
      childMeshes[i].flattenTree(_workingT,L);
  }//endfunction
  
  /**
   * renders the given mesh tree
   */
  static void render(RenderingContext gl,Mesh tree,int viewWidth, int viewHeight, double aspect) 
  {
    gl.clearColor(0.0, 0.0, 0.0, 1.0);  // clear color and alpha
    
    // Basic viewport setup and clearing of the screen
    gl.viewport(0, 0, viewWidth, viewHeight);
    gl.clear(COLOR_BUFFER_BIT | DEPTH_BUFFER_BIT);
    gl.enable(DEPTH_TEST);
    gl.enable(BLEND);
    
    // ----- calculate lighting info 
    Float32List lPs = new Float32List(lightPts.length*3);     // light positions
    Float32List lCs = new Float32List(lightPts.length*3);     // light colors
    for (int i=lightPts.length-1; i>-1; i--)
    {
      PointLight lpt = lightPts[i];
      lPs[i*3+0] = lpt.px*viewT.aa+lpt.py*viewT.ab+lpt.pz*viewT.ac+viewT.ad;  // transformed lightPt x
      lPs[i*3+1] = lpt.px*viewT.ba+lpt.py*viewT.bb+lpt.pz*viewT.bc+viewT.bd;  // transformed lightPt y
      lPs[i*3+2] = lpt.px*viewT.ca+lpt.py*viewT.cb+lpt.pz*viewT.cc+viewT.cd;  // transformed lightPt z
      lCs[i*3+0] = lpt.r;
      lCs[i*3+1] = lpt.g;
      lCs[i*3+2] = lpt.b;
    }//endfor
    
    List<Mesh> Mshs = new List<Mesh>();
    tree.flattenTree(viewT,Mshs);
    
    for (int i=Mshs.length-1; i>-1; i--)
    {
      Mesh m=Mshs[i];   // current mesh to render
      if (m.idxBuffer!=null && m.vertBuffer!=null)
      {
        if (m.prepOnState!=stateId) m.prepareForRender();         // recompile shader code if state changed
        
        gl.useProgram(m.prog);                                    // set shader program
        gl.bindBuffer(ELEMENT_ARRAY_BUFFER, m.idxBuffer);         // bind index buffer
        gl.bindBuffer(ARRAY_BUFFER, m.vertBuffer);                // bind vertex buffer
        int aLoc = gl.getAttribLocation(m.prog, 'aPosn');
        gl.enableVertexAttribArray(aLoc);
        gl.vertexAttribPointer(aLoc, 3, FLOAT, false, 8*4, 0);    // specify vertex
        aLoc = gl.getAttribLocation(m.prog, 'aUV');
        gl.enableVertexAttribArray(aLoc);
        gl.vertexAttribPointer(aLoc, 2, FLOAT, false, 8*4, 3*4);  // specify UV
        aLoc = gl.getAttribLocation(m.prog, 'aNorm');
        gl.enableVertexAttribArray(aLoc);
        gl.vertexAttribPointer(aLoc, 3, FLOAT, false, 8*4, 5*4);  // specify normal
        
        gl.uniformMatrix4fv(gl.getUniformLocation(m.prog, 'uPMatrix'), false, persT.buf());
        gl.uniformMatrix4fv(gl.getUniformLocation(m.prog, 'uMVMatrix'), false, m._workingT.buf());
        
        gl.uniform3fv(gl.getUniformLocation(m.prog, 'lPoints'),lPs);  // upload light posns
        gl.uniform3fv(gl.getUniformLocation(m.prog, 'lColors'),lCs);  // upload light colors
        
        gl.drawArrays(TRIANGLES, 0, m.numTris*3);
        //gl.drawElements(TRIANGLES, m.numTris*3, UNSIGNED_SHORT, 0);
      }//endif
    }//endfor
    
  }//endfunction
  
  /**
   * compiles and uploads shader program
   */
  static Program uploadShaderProgram(String vertSrc, String fragSrc) 
  {
    Shader vertShader = gl.createShader(VERTEX_SHADER);
    gl.shaderSource(vertShader, vertSrc);
    gl.compileShader(vertShader);

    Shader fragShader = gl.createShader(FRAGMENT_SHADER);
    gl.shaderSource(fragShader, fragSrc);
    gl.compileShader(fragShader);

    Program prog = gl.createProgram();
    gl.attachShader(prog, vertShader);
    gl.attachShader(prog, fragShader);
    gl.linkProgram(prog);

    Object linkStat = gl.getProgramParameter(prog, LINK_STATUS);
    print("linkStatus=" + linkStat.toString());

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
    persT = Matrix4.perspective(fov, aspectRatio, nearClip, farClip);
    viewT = getViewTransform(px,py,pz,tx,ty,tz);
    camT = viewT.inverse();
    return camT.scale(1.0,1.0,1.0); // duplicate and return
  }//endfunction
  
  /**
   * sets point lights positions and colors
   */
  static void setLights(List<PointLight> lights)
  {
    if (lightPts.length!=lights.length) stateId++; // to trigger recompile of shader codes
    lightPts=lights;
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
    return new Mesh(VData);
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
                    0,0,0,  // normal a
                    V[I[i+1]*3+0],V[I[i+1]*3+1],V[I[i+1]*3+2],  // vertex b
                    U[(i*2)%ul+2],U[(i*2)%ul+3],
                    0,0,0,  // normal b
                    V[I[i+2]*3+0],V[I[i+2]*3+1],V[I[i+2]*3+2],  // vertex c
                    U[(i*2)%ul+4],U[(i*2)%ul+5],
                    0,0,0]); // normal c
    }
    return new Mesh(VData);
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
    //M.compressGeometry();
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
    //M.compressGeometry();
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
    HttpRequest.getString(url).then((String s) {handle(Mesh.parseObj(s));});
  }//endfunction
    
  /**
   * parses a given obj format string data s to mesh with null texture
   * Mtls: [id1,bmd1,id2,bmd2,...]
   */
  static Mesh parseObj(String s)
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
        mmesh.addChild(cm);
      }
    }//endfor g
    
    return mmesh; // returns mesh with submeshes in it
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