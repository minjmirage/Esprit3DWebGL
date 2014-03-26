import 'dart:typed_data';
import 'dart:math';

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

