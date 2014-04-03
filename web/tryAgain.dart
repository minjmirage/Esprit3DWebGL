import 'dart:html';
import 'dart:web_gl';
import 'dart:math';
import 'core3d/mesh.dart';

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

  loadImgs(["weaveDiff.jpg","weaveNorm.jpg","weaveSpec.jpg","hexDiff.jpg","hexNorm.jpg","hexSpec.jpg"],(List<ImageElement> Texs)
      {
        scene.setTexture(Texs[0],true);
        scene.setNormalMap(Texs[1],true);
        scene.setSpecularMap(Texs[2],true);
        Mesh cube = Mesh.createCube(-13.0,-13.0,-13.0,true);
        cube.setTexture(Texs[3]);
        cube.setNormalMap(Texs[4]);
        cube.setSpecularMap(Texs[5]);
        scene.addChild(cube);
        tick(0.0);  // start main loop!
      });
}//endmain

/**
 * convenience function to load a list of images
 */
void loadImgs(List<String> URLs,callBack(List<ImageElement> Texs),[List<ImageElement> Texs=null])
{
  if (Texs==null) Texs=new List<ImageElement>();
  print("Loading "+Texs.length.toString()+" : "+URLs[Texs.length]);
  ImageElement img = new ImageElement(src:URLs[Texs.length]);
  img.onLoad.listen((dat)
      {
        Texs.add(img);
        if (Texs.length<URLs.length)
          loadImgs(URLs,callBack,Texs);
        else
          callBack(Texs);
      });
}//endfunction

tick(double time)
{
  window.animationFrame.then(tick);
  print("time:"+time.toString());
  double r = 5+sin(time/500);
  double sinT = sin(sin(time/1000));

  Mesh.setLights([new PointLight(10.0*sin(time/500),0.0,10.0*cos(time/500),1.0,1.0,1.0)]);
                 // new PointLight(10.0*sin(time/500+PI*2/3),0.0,10.0*cos(time/500+PI*2/3),1.0,0.0,1.0),
                 // new PointLight(10.0*sin(time/500+PI*4/3),0.0,10.0*cos(time/500+PI*4/3),0.0,1.0,1.0)]);
  Mesh.setCamera(r*cos(sinT/10)*sin(time/5000), r*sin(sinT/10), r*cos(sinT/10)*cos(time/5000), 0.0,0.0,0.0);

  Mesh.render(gl,scene,800, 600, 800 / 600);
}

