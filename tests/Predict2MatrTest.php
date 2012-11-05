<?php defined('SYSPATH') OR die('Kohana bootstrap needs to be included before tests run');

class Predict2MatrTest extends Unittest_TestCase{
  function providerForJob021(){
    return array(
      array(//matrix  around X axis
        'n'=>1,'a'=>deg2rad(36.0),
        'p'=>Array(null,
          Array(null,1.0, 0.0,0.0),
          Array(null,0.0, 0.809016994374947,0.587785252292473),
          array(null,0.0,-0.587785252292473,0.809016994374947)
        ),
      ),array(//matrix  around Y axis
        'n'=>2,'a'=>deg2rad(36.0),
        'p'=>Array(null,
          array(null,0.809016994374947,0.0, -0.587785252292473),
          array(null,0.0,1.0,0.0),
          array(null,0.587785252292473,0.0,0.809016994374947),
        ),
      ),array(//matrix  around Z axis
        'n'=>3,'a'=>deg2rad(36.0),
        'p'=>Array(null,
          array(null,0.809016994374947,0.587785252292473,0.0),
          array(null,-0.587785252292473,0.809016994374947,0.0),
          array(null,0.0,0.0,1.0)    
        ),
      )
    );
  }
  /**
  * 2.1 To calculate matrixes of rotation
  * @dataProvider providerForJob021
  */
  function testForJob021($n,$a,$p1){
    $p=Predict2_Matr::ForRotMatr($n,$a);
    $p[0]=$p[1][0]=$p[2][0]=$p[3][0]=null;//wrap dimensions
    $this->assertEquals($p,$p1,'',1e-7);
  }
  function providerForJob022(){
    return array(
      array(
      'matrix'=>array(null,//output
        array(null,-0.186639361025995,  -0.960176274304416,  -0.207911690817759),
        array(null,-0.888091727606158,   0.255385612942791,  -0.382192715863795),
        array(null, 0.420070032581851,   0.113312448410253,  -0.900389613868328),
      ),
      'vector_r'=>array(//input
        1=>381794.865023,2=>136522.450954,3=>24396.993002
      ),
      'vector_v'=>array(//output
        1=>-207415.988033,2=>-313527.344466,3=>153883.475461
      )
    ));
  }
  
  /**
  * 2.2 To multiply matrix by vector
  * @dataProvider providerForJob022
  */
  function testForJob022($p1,$r,$v1){
    $a=Predict2_Matr::ForRotMatr(1,deg2rad(23.0)); 
    $b=Predict2_Matr::ForRotMatr(2,deg2rad(168.0)); 
    $c=Predict2_Matr::ToMultMatr($a,$b); 
    $a=Predict2_Matr::ForRotMatr(3,deg2rad(79.0)); 
    $p=Predict2_Matr::ToMultMatr($c,$a);
    $p[0]=$p[1][0]=$p[2][0]=$p[3][0]=null;//wrap dimensions
    $this->assertEquals($p,$p1,'',1e-7);
   
    $v=Predict2_Matr::MultMatrVec($p,$r);
    $this->assertEquals($v,$v1,'',1e-6);
  }  

  function providerForJob023(){
    return array(
      array(//input
      'matrix A'=>array(
      1=>array(null,0.999987104372264,  -0.004657817791984,  -0.002023813872801),
      2=>array(null,0.004657817791642,   0.999989152296765,  -0.000004713477259),
      3=>array(null,0.002023813873587,  -0.000004713139788,   0.999997952075499),
      ),
      'matrix B'=>array(
      1=>array(null, 0.999999995979347,   0.000082275361072,   0.000035666114213),
      2=>array(null,-0.000082275089649,   0.999999996586437,  -0.000007611504096),
      3=>array(null,-0.000035666740330,   0.000007608569633,   0.999999999334997),
      ),
      'matrix C=A*B'=>array(
      1=>array(null,0.999987555756883,  -0.004575558874333,  -0.001988112764180),
      2=>array(null,0.004575543743875,   0.999989532071017,  -0.000012158772523),
      3=>array(null,0.001988147585937,   0.000003061924296,   0.999998023627948),
      ))
    );
  }
  
  /**
  * 2.3 To multiply matrix by matrix
  * @dataProvider providerForJob023
  */
  function testForJob023($a,$b,$c1){
    $c=Predict2_Matr::ToMultMatr($a,$b);
    $c[1][0] = $c[2][0]=$c[3][0]=null;//wrap
    $this->assertEquals($c,$c1,'',1e-8);
  }

  function providerForJob024(){
    return array(
      array(//input
      'matrix A'=>array(
      1=>array(null,0.999987555756883, -0.004575558874333, -0.001988112764180),
      2=>array(null,0.004575543743875,  0.999989532071017, -0.000012158772523),
      3=>array(null,0.001988147585937,  0.000003061924296,  0.999998023627948),
      ),
      'matrix B'=>array(
      1=>array(null, 0.999987555756883,  0.004575543743875,  0.001988147585937),
      2=>array(null,-0.004575558874333,  0.999989532071017,  0.000003061924296),
      3=>array(null,-0.001988112764180, -0.000012158772523,  0.999998023627948),
      ),
      'matrix C=A*B'=>array(
      1=>array(null, 1.000000000000000, -0.000000000000000,  0.000000000000000),
      2=>array(null,-0.000000000000000,  1.000000000000000,  0.000000000000000),
      3=>array(null, 0.000000000000000,  0.000000000000000,  1.000000000000000),
      ))
    );
  }
  
  /**
  * 2.4 To transponir matrix
  * @dataProvider providerForJob024
  */
  function testForJob024($a,$b1,$c1){
    $b=Predict2_Matr::TranspMatr($a);
    $c=Predict2_Matr::ToMultMatr($a,$b);
    $b[1][0]=$b[2][0]=$b[3][0]=null;
    $c[1][0]=$c[2][0]=$c[3][0]=null;
    $this->assertEquals($b,$b1,'',1e-8);
    $this->assertEquals($c,$c1,'',1e-8);
  }

  function providerForJob025(){
    return array(
      array(//input 
      'spheric position'=>array(
      'l'=> 3.787930671498,//   longitude in radian
      'b'=>-0.000003961191,//   latitude in radian
      'r'=> 148552579.338,//     range in km
      ),
      'rectangular position'=>array(
      1=> -118588727.486,//  x position in km
      2=>  -89468332.616,//  y position in km
      3=>       -588.445//  z position in km
      ))
    );
  }
  
  /**
  * 2.5 From spheric to rectangular position
  * @dataProvider providerForJob025
  */
  function testForJob025($s,$v1){
    $v=Predict2_Matr::DescFromSpherCoor($s['l'],$s['b'],$s['r']); 
    $this->assertEquals($v,$v1,'',1e-3);
  }
  
}