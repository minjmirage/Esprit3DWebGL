import 'dart:html';
import 'dart:web_gl';
import 'dart:typed_data';

import 'matrix4.dart';

CanvasElement canvas = querySelector("#glCanvas");
RenderingContext gl;
Matrix4 mvMatrix, pMatrix;

void main() 
{
  // Nab the context we'll be drawing to.
  print("canvas:" + canvas.width.toString() + "x" + canvas.height.toString());
  gl = canvas.getContext3d();
  if (gl == null) return;

  Mesh cube = Mesh.createCube(0.1,0.2,0.3);
  Mesh.render(gl,cube,800, 600, 800 / 600);
}

/**
 * Statically draw a triangle and a square!
 */
class Mesh 
{
  Program prog;           // the render program for this mesh
  Buffer vertBuffer;      // uploaded vertices data to GPU
  Buffer idxBuffer;       // uploaded indices data to GPU
  List<Mesh> childMeshes; // list of children
  Matrix4 transform;      // transform for this mesh
  int numTris;            // num of triangles to render

  /**
   * creates new mesh from given geometry data
   */
  Mesh([List<double> vertData=null,List<int> idxData=null])
  {
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
    
    // ----- create and upload program
    String vertSrc = '''
        attribute vec3 aPosn;
        attribute vec2 aUV;
        attribute vec3 aNorm;
        
        varying vec3 vPosn;
        varying vec2 vUV;
        varying vec3 vNorm;
        
        uniform mat4 uMVMatrix;
        uniform mat4 uPMatrix;
        uniform mat3 uNMatrix;
    
        void main(void) {
            vPosn = (uMVMatrix * vec4(aPosn, 1.0)).xyz;
            gl_Position = uPMatrix * vec4(vPosn, 1.0);
            vUV = aUV;
            vNorm = (uMVMatrix * vec4(aNorm, 0.0)).xyz;
        }
      ''';
    String fragSrc = '''
        precision highp float;
        
        varying vec2 vUV;
        varying vec3 vNorm;
        varying vec3 vPosn;
               
        void main(void) {
            gl_FragColor = vec4(1.0,1.0,1.0,1.0);
        }
      '''; 
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
    prog = uploadShaderProgram(vertSrc,fragSrc);
    gl.useProgram(prog);
  }//endconstr
  
  /**
   * adds given mesh to this render tree
   */
  void addChild(Mesh m)
  {
    childMeshes.add(m);
  }//endfunction
  
  /**
   * renders the given mesh tree
   */
  static void render(RenderingContext gl,Mesh m,num viewWidth, num viewHeight, num aspect) 
  {
    gl.clearColor(0.0, 0.0, 0.0, 1.0);  // clear color and alpha
    
    // Basic viewport setup and clearing of the screen
    gl.viewport(0, 0, viewWidth, viewHeight);
    gl.clear(COLOR_BUFFER_BIT | DEPTH_BUFFER_BIT);
    gl.enable(DEPTH_TEST);
    gl.enable(BLEND);
    
    pMatrix = Matrix4.perspective(45.0, aspect, 0.1, 100.0);
    mvMatrix = new Matrix4().rotX(1.0).translate(0.0, 0.0, -7.0);
        
    print("pmatrix="+pMatrix.toString());
    print("mvmatrix="+mvMatrix.buf().toString());
        
    // Here's that bindBuffer() again, as seen in the constructor
    gl.bindBuffer(ELEMENT_ARRAY_BUFFER, m.idxBuffer);   // bind index buffer
    gl.bindBuffer(ARRAY_BUFFER, m.vertBuffer);          // bind vertex buffer
    int aLoc = gl.getAttribLocation(m.prog, 'aPosn');
    gl.enableVertexAttribArray(aLoc);
    gl.vertexAttribPointer(aLoc, 3, FLOAT, false, 0, 0);
    aLoc = gl.getAttribLocation(m.prog, 'aNorm');
    gl.enableVertexAttribArray(aLoc);
    gl.vertexAttribPointer(aLoc, 3, FLOAT, false, 0, 3*4);
    aLoc = gl.getAttribLocation(m.prog, 'aUV');
    gl.enableVertexAttribArray(aLoc);
    gl.vertexAttribPointer(aLoc, 3, FLOAT, false, 0, 6*4);
    gl.uniformMatrix4fv(gl.getUniformLocation(m.prog, 'uPMatrix'), false, pMatrix.buf());
    gl.uniformMatrix4fv(gl.getUniformLocation(m.prog, 'uMVMatrix'), false, mvMatrix.buf());
    gl.drawElements(TRIANGLES, m.numTris*3, UNSIGNED_SHORT, 0);
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
                   3,7,2, 2,7,6];// bottom
                       
    List<double> U = [0.0,1.0, 0.0,0.0, 1.0,1.0, 1.0,1.0, 0.0,0.0, 1.0,0.0];
    
    int i=0;
    int ul=U.length;
    List<num> VData = new List<double>();
    if (soft)
    {
      for (i=0; i<I.length; i+=3)
      VData.addAll([V[I[i+0]*3+0],V[I[i+0]*3+1],V[I[i+0]*3+2],  // vertex a
                    V[I[i+0]*3+0],V[I[i+0]*3+1],V[I[i+0]*3+2],  // normal a
                    U[(i*2)%ul+0],U[(i*2)%ul+1],
                    V[I[i+1]*3+0],V[I[i+1]*3+1],V[I[i+1]*3+2],  // vertex b
                    V[I[i+1]*3+0],V[I[i+1]*3+1],V[I[i+1]*3+2],  // normal b
                    U[(i*2)%ul+2],U[(i*2)%ul+3],
                    V[I[i+2]*3+0],V[I[i+2]*3+1],V[I[i+2]*3+2],  // vertex c
                    V[I[i+2]*3+0],V[I[i+2]*3+1],V[I[i+2]*3+2],  // normal c
                    U[(i*2)%ul+4],U[(i*2)%ul+5]]);
    }
    else
    {
      for (i=0; i<I.length; i+=3)
      VData.addAll([V[I[i+0]*3+0],V[I[i+0]*3+1],V[I[i+0]*3+2],  // vertex a
                    0,0,0,  // normal a
                    U[(i*2)%ul+0],U[(i*2)%ul+1],
                    V[I[i+1]*3+0],V[I[i+1]*3+1],V[I[i+1]*3+2],  // vertex b
                    0,0,0,  // normal b
                    U[(i*2)%ul+2],U[(i*2)%ul+3],
                    V[I[i+2]*3+0],V[I[i+2]*3+1],V[I[i+2]*3+2],  // vertex c
                    0,0,0,  // normal c
                    U[(i*2)%ul+4],U[(i*2)%ul+5]]);
    }
    return new Mesh(VData);
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
    List<num> V = new List<num>();    // list to contain vertices data
    List<num> T = new List<num>();    // list to contain texture coordinates data
    List<num> N = new List<num>();    // list to contain normals data
    
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
          V.add(num.parse(A[j]));
      }
      else if (D[i].startsWith('vt '))   // ----- if vertex uv definition
      {
        A = (D[i].substring(2)).split(' ');
        for (j=A.length-1; j>=0; j--)
          if (A[j]=="") A.removeAt(j);
        for (j=0; j<A.length && j<2; j++)   // restrict to u,v instead of u,v,t
          T.add(num.parse(A[j]));
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
          N.add(num.parse(A[j]));
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
    
    num onNumParseError(String s) {return 1;};
    
    for (int g=0; g<G.length; g++)
    {
      F = G[g];
      
      // ----- import faces data -----------------------------
      List<num> verticesData = new List<num>();  // to contain [vx,vy,vz,nx,ny,nz,u,v, ....]
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
            List p = f[j];      // data of a point: [v,uv,n] to num
            for (int k=0; k<p.length; k++)
              p[k] = (num.parse(p[k],onNumParseError)).toInt()-1; 
          }
          
          // ----- triangulate higher order polygons
          while (f.length>=3)
          {
            A = new List();
            for (j=0; j<3; j++)
              A.addAll(f[j]);
            // A: [v,uv,n,v,uv,n,v,uv,n]
            
            // ----- get vertices --------------------------------
            num vax = V[A[0]*3+0];
            num vay = V[A[0]*3+1];
            num vaz = V[A[0]*3+2];
            num vbx = V[A[3]*3+0];
            num vby = V[A[3]*3+1];
            num vbz = V[A[3]*3+2];
            num vcx = V[A[6]*3+0];
            num vcy = V[A[6]*3+1];
            num vcz = V[A[6]*3+2];
            
            // ----- get normals ---------------------------------
            num px = vbx - vax;
            num py = vby - vay;
            num pz = vbz - vaz;
              
            num qx = vcx - vax;
            num qy = vcy - vay;
            num qz = vcz - vaz;
            // normal by determinant
            num nx = py*qz-pz*qy;  //  unit normal x for the triangle
            num ny = pz*qx-px*qz;  //  unit normal y for the triangle
            num nz = px*qy-py*qx;  //  unit normal z for the triangle
            
            num nax = nx;   // calculated normals
            num nay = ny;
            num naz = nz;
            num nbx = nx;
            num nby = ny;
            num nbz = nz;
            num ncx = nx;
            num ncy = ny;
            num ncz = nz;
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
            num ua = 0;
            num va = 0;
            num ub = 1;
            num vb = 0;
            num uc = 0;
            num vc = 1;
            if (T.length>0)
            {
              ua = T[A[1]*2+0];
              va = 1-T[A[1]*2+1];
              ub = T[A[4]*2+0];
              vb = 1-T[A[4]*2+1];
              uc = T[A[7]*2+0];
              vc = 1-T[A[7]*2+1];
            }
                                  
            verticesData.addAll([ vax,vay,vaz, nax,nay,naz, ua,va,
                                  vbx,vby,vbz, nbx,nby,nbz, ub,vb,
                                  vcx,vcy,vcz, ncx,ncy,ncz, uc,vc]);
                      
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