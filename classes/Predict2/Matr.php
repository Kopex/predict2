<?php defined('SYSPATH') OR die('Kohana bootstrap needs to be included before tests run');
/**
 * @package     Predict2
 * @copyright   (c) 2012 Kopex
 */

class Predict2_Matr{

  /**
  * Convert position velosity from FixEqu to TruEqu
  * 
  * @param float $tb MJD
  * @param Predict2_Model_Vector $rb x,y,z in km
  * @param Predict2_Model_Vector $vb vx,vy,vz in km/s
  * @param Predict2_Model_Vector $ra x,y,z in km
  * @param Predict2_Model_Vector $va vx,vy,vz in km/s
  */
  static function PosVelToTruEqu($tb,Predict2_Model_Vector $rb,Predict2_Model_Vector $vb,
    Predict2_Model_Vector $ra,Predict2_Model_Vector $va){
    self::FromFixToTrueEquM($tb,$rm); // simple matrix UnForRom 
    self::MultMatrVec(rm,rb,ra);    // to true equator system UnForFun 
    self::MultMatrVec(rm,vb,va);    // to true equator system 
  }
  /**
  * Convert position velosity from TruEqu to FixEqu
  * 
  * @param mixed $tb
  * @param Predict2_Model_Vector $ra
  * @param Predict2_Model_Vector $va
  * @param Predict2_Model_Vector $rb
  * @param Predict2_Model_Vector $vb
  */
  static function PosVelToFixEqu($tb,Predict2_Model_Vector $ra,Predict2_Model_Vector $va,
    Predict2_Model_Vector $rb,Predict2_Model_Vector $vb){
    self::FromTrueToFixEquM($tb,$rm); // simple matrix UnForRom 
    self::MultMatrVec($rm,$ra,$rb);    // to fixed equator system UnForFun 
    self::MultMatrVec($rm,$va,$vb);    // to fixed equator system 
  }
  static function FromTrueToFixEquM ($tc,$pn){
    self::FromFixToTrueEquM($tc,$pun); // matrix from fixed to true equator 
    self::TranspMatr(pun,pn);        // to transponir UnForFun 
  }
    
  static function FromFixToTrueEquM ( $tc,$pn){
    self::ClcPrecMatr(MJD2000,$tc,$prec); // UnForPnm precession matrix
    self::ClcNut($tc,$dpsi,$deps) ;  // obliquity radian UnForNut }
    self::ClcNutMatr(tc,dpsi,deps,ToGetEpsMean(tc),nutm); // UnForPnm UnForNut 
    self::ToMultMatr(nutm,prec,pn); // from fixed to true equator 
  }

  /* to multiple to simple 3*3 matrix */
  static function ToMultMatr($ma,$mb){
    for (  $i = 1 ;  $i <= 3 ; $i++ )
     for (  $j = 1 ;  $j <= 3 ; $j++ ){
         $mc[$i][$j] = 0.0;
         for (  $k = 1 ;  $k <= 3 ; $k++ )
           $mc[$i][$j] = $mc[$i][$j]+$ma[$i][$k]*$mb[$k][$j];
     }
    return $mc;
  }

  static function MultMatrVec( $p,$v){
    For($i=1;$i<=3;$i++){
       $vp[$i]=0.0;
       For($j=1;$j<=3;$j++)
         $vp[$i]=$vp[$i]+$p[$i][$j]*$v[$j];
    }
    return $vp;
  }
  static function TranspMatr($ma){
    For($i=1;$i<=3;$i++){
      For($j=1;$j<=3;$j++)
        $mb[$i][$j]=$ma[$j][$i];
    }
  return $mb;
  }
/* ClcPrecMatr is a procedure
    to calculate the Matrix of precession
       for change from   Epoch  to  Tmjd  or another words
       from the fixed equator and ecliptic to the instantaneous ones.
       Fundamental  EPOCH  is  J2000.0
         but  Epoch  as parameter
              may be any fixed epoch in  MJD
                 and  Tmjd  is epoch of the Date in  MJD .
       PrecMatr[1..3,1..3] is the resulting massive .
       global const - SecRad - for change from arcsec to radians . */
  static function ClcPrecMatr ( $Epoch , $Tmjd, $PrecMatr){
    self::ClcPrecAngles($Epoch,$Tmjd,$Dzita0,$Teta,$Zet); // angles in ArcSec 
    self::ForRotMatr(3,-SecRad*Dzita0,aa); // procedure from UnForFun 
    self::ForRotMatr(2,SecRad*Teta,ab);    // procedure for rotation 
    self::ForRotMatr(3,-SecRad*Zet,ac);
    self::ToMultMatr(ab,aa,am); // procedure from UnForFun 
    self::ToMultMatr(ac,am,PrecMatr); // this is our result 
  }
 
  /* fill the elements of the simple rotation matrix */
  static function ForRotMatr ( $iaxis , $angle){
    $c = cos($angle);
    $s = sin($angle);
    for (  $i = 1 ;  $i <= 3 ; $i++ )
     for (  $j = 1 ;  $j <= 3 ; $j++ )
       $matr[$i][$j] = ( $i == $j )?1.0:0.0; /* unit matrix */
    switch ( $iaxis ){
     case 1 : /* rotate around ox axis */
           $matr[2][2] = +$c;
           $matr[2][3] = +$s;
           $matr[3][2] = -$s;
           $matr[3][3] = +$c;
           break;
     case 2 : /* rotate around oy axis */
           $matr[1][1] = +$c;
           $matr[1][3] = -$s;
           $matr[3][1] = +$s;
           $matr[3][3] = +$c;
           break;
     case 3 :
           $matr[1][1] = +$c;
           $matr[1][2] = +$s;
           $matr[2][1] = -$s;
           $matr[2][2] = +$c;
           break;
    } /* case for axis of rotation */
    return $matr;
  }

  static function ClcNutMatr ( $Tmjd,$DeltaPsi,$DeltaEps,$EpsMean,$MatrNut){
    self::ForRotMatr(1,EpsMean,MatrEps);   // procedure from UnForFun 
    self::ForRotMatr(3,-DeltaPsi,MatrPsi); // matrix of rotation 
    self::ToMultMatr(MatrPsi,MatrEps,MatrCur); // procedure from UnForFun 
    self::ForRotMatr(1,-EpsMean-DeltaEps,MatrEps);
    self::ToMultMatr(MatrEps,MatrCur,MatrNut); // this is our result 
  }
  /* clcnut is a void
     to calculate the parameters of nutation in longitude and obliquity
        and value of the mean obliquity in radians on the moment tmjd.
        the arguments and coefficients
            of nutation theory are contained in the nutarg and coefnut .
        global parameter - secrad - for change from arcsecs to radians . */
   static function ClcNut ($tmjd , &$deltapsi , &$deltaeps ){
     /* 1980 iau theory of nutation  ( J. M. Wahr ) */
     //[106][5] /* arguments */
     $nutarg = array(   
       array(1 ,0 ,0 ,0 ,0),  array(2 ,0 ,0 ,0 ,0), array(1 ,-2 ,0 ,2 ,0), array(0 ,2 ,0 ,-2 ,0), array(2 ,-2 ,0 ,2 ,0),array(0 ,1 ,-1 ,0 ,-1),
       array(1 ,0 ,-2 ,2 ,-2),array(1 ,2 ,0 ,-2 ,0),array(2 ,0 ,0 ,2 ,-2), array(0 ,0 ,1 ,0 ,0),  array(2 ,0 ,1 ,2 ,-2),array(2 ,0 ,-1 ,2 ,-2),
       array(1 ,0 ,0 ,2 ,-2), array(0 ,2 ,0 ,0 ,-2),array(0 ,0 ,0 ,2 ,-2), array(0 ,0 ,2 ,0 ,0),  array(1 ,0 ,1 ,0 ,0), array(2 ,0 ,2 ,2 ,-2),
       array(1 ,0 ,-1 ,0 ,0), array(1 ,-2 ,0 ,0 ,2),array(1 ,0 ,-1 ,2 ,-2),array(1 ,2 ,0 ,0 ,-2), array(1 ,0 ,1 ,2 ,-2),array(0 ,1 ,0 ,0 ,-1),
       array(0 ,2 ,1 ,0 ,-2), array(1 ,0 ,0 ,-2 ,2),array(0 ,0 ,1 ,-2 ,2), array(2 ,0 ,1 ,0 ,0),  array(1 ,-1 ,0 ,0 ,1),array(0 ,0 ,1 ,2 ,-2),
       array(2 ,0 ,0 ,2 ,0),  array(0 ,1 ,0 ,0 ,0), array(1 ,0 ,0 ,2 ,0),  array(2 ,1 ,0 ,2 ,0),  array(0 ,1 ,0 ,0 ,-2),array(2 ,-1 ,0 ,2 ,0),
       array(0 ,0 ,0 ,0 ,2),  array(1 ,1 ,0 ,0 ,0), array(1 ,-1 ,0 ,0 ,0), array(2 ,-1 ,0 ,2 ,2), array(1 ,1 ,0 ,2 ,0), array(2 ,0 ,0 ,2 ,2),
       array(0 ,2 ,0 ,0 ,0),  array(2 ,1 ,0 ,2 ,-2),array(2 ,2 ,0 ,2 ,0),  array(0 ,0 ,0 ,2 ,0),  array(1 ,-1 ,0 ,2 ,0),array(1 ,-1 ,0 ,0 ,2),
       array(1 ,1 ,0 ,0 ,-2), array(1 ,-1 ,0 ,2 ,2),array(0 ,1 ,1 ,0 ,-2), array(2 ,0 ,1 ,2 ,0),  array(2 ,0 ,-1 ,2 ,0),array(2 ,1 ,0 ,2 ,2),
       array(0 ,1 ,0 ,0 ,2),  array(2 ,2 ,0 ,2 ,-2),array(1 ,0 ,0 ,0 ,2),  array(1 ,0 ,0 ,2 ,2),  array(1 ,1 ,0 ,2 ,-2),array(1 ,0 ,0 ,0 ,-2),
       array(0 ,1 ,-1 ,0 ,0), array(1 ,2 ,0 ,2 ,0), array(0 ,0 ,1 ,0 ,-2), array(0 ,1 ,0 ,-2 ,0), array(0 ,0 ,0 ,0 ,1), array(0 ,1 ,1 ,0 ,0),
       array(0 ,1 ,0 ,2 ,0),  array(2 ,1 ,-1 ,2 ,0),array(2 ,-1 ,-1 ,2 ,2),array(1 ,-2 ,0 ,0 ,0), array(2 ,3 ,0 ,2 ,0), array(2 ,0 ,-1 ,2 ,2),
       array(2 ,1 ,1 ,2 ,0),  array(1 ,-1 ,0 ,2 ,-2),array(1 ,2 ,0 ,0 ,0), array(2 ,1 ,0 ,0 ,0),  array(0 ,3 ,0 ,0 ,0), array(2 ,0 ,0 ,2 ,1),
       array(2 ,-1 ,0 ,0 ,0), array(0 ,1 ,0 ,0 ,-4),array(2 ,-2 ,0 ,2 ,2), array(2 ,-1 ,0 ,2 ,4), array(0 ,2 ,0 ,0 ,-4),array(2 ,1 ,1 ,2 ,-2),
       array(1 ,1 ,0 ,2 ,2),  array(2 ,-2 ,0 ,2 ,4),array(2 ,-1 ,0 ,4 ,0), array(0 ,1 ,-1 ,0 ,-2),array(1 ,2 ,0 ,2 ,-2),array(2 ,2 ,0 ,2 ,2),
       array(1 ,1 ,0 ,0 ,2),  array(2 ,0 ,0 ,4 ,-2),array(2 ,3 ,0 ,2 ,-2), array(0 ,1 ,0 ,2 ,-2), array(1 ,0 ,1 ,2 ,0), array(1 ,-1 ,-1 ,0 ,2),
       array(1 ,0 ,0 ,-2 ,0), array(2 ,0 ,0 ,2 ,-1),array(0 ,0 ,1 ,0 ,2),  array(0 ,1 ,0 ,-2 ,-2),array(1 ,0 ,-1 ,2 ,0),array(1 ,1 ,1 ,0 ,-2),
       array(0 ,1 ,0 ,-2 ,2), array(0 ,2 ,0 ,0 ,2), array(2 ,0 ,0 ,2 ,4),  array(0 ,0 ,1 ,0 ,1)
     );
     /* 1980 iau theory of nutation  ( J. M. Wahr ) 
     [106][4]amplitudes */
     $coefnut =array(  
        array( -17.199600 ,-0.017420 ,9.202500 ,0.000890), array( 0.206200 ,0.000020 ,-0.089500 ,0.000050),array( 0.004600 ,0.000000 ,-0.002400 ,0.000000),   array( 0.001100 ,0.000000 ,0.000000 ,0.000000),
        array(-0.000300 ,0.000000 ,0.000100 ,0.000000),    array(-0.000300 ,0.000000 ,0.000000 ,0.000000), array(-0.000200 ,0.000000 ,0.000100 ,0.000000),    array( 0.000100 ,0.000000 ,0.000000 ,0.000000),
        array(-1.318700 ,-0.000160 ,0.573600 ,-0.000310),  array( 0.142600 ,-0.000340 ,0.005400 ,-0.000010),array(-0.051700 ,0.000120 ,0.022400 ,-0.000060),   array( 0.021700 ,-0.000050 ,-0.009500 ,0.000030),
        array( 0.012900 ,0.000010 ,-0.007000 ,0.000000),   array( 0.004800 ,0.000000 ,0.000100 ,0.000000),  array(-0.002200 ,0.000000 ,0.000000 ,0.000000),    array( 0.001700 ,-0.000010 ,0.000000 ,0.000000),
        array(-0.001500 ,0.000000 ,0.000900 ,0.000000),    array(-0.001600 ,0.000010 ,0.000700 ,0.000000),  array(-0.001200 ,0.000000 ,0.000600 ,0.000000),    array(-0.000600 ,0.000000 ,0.000300 ,0.000000),
        array(-0.000500 ,0.000000 ,0.000300 ,0.000000),    array( 0.000400 ,0.000000 ,-0.000200 ,0.000000), array( 0.000400 ,0.000000 ,-0.000200 ,0.000000),   array(-0.000400 ,0.000000 ,0.000000 ,0.000000),
        array( 0.000100 ,0.000000 ,0.000000 ,0.000000),    array( 0.000100 ,0.000000 ,0.000000 ,0.000000),  array(-0.000100 ,0.000000 ,0.000000 ,0.000000),    array( 0.000100 ,0.000000 ,0.000000 ,0.000000),
        array( 0.000100 ,0.000000 ,0.000000 ,0.000000),    array(-0.000100 ,0.000000 ,0.000000 ,0.000000),  array(-0.227400 ,-0.000020 ,0.097700 ,-0.000050),  array( 0.071200 ,0.000010 ,-0.000700 ,0.000000),
        array(-0.038600 ,-0.000040 ,0.020000 ,0.000000),   array(-0.030100 ,0.000000 ,0.012900 ,-0.000010), array(-0.015800 ,0.000000 ,-0.000100 ,0.000000),   array( 0.012300 ,0.000000 ,-0.005300 ,0.000000),
        array( 0.006300 ,0.000000 ,-0.000200 ,0.000000),   array( 0.006300 ,0.000010 ,-0.003300 ,0.000000), array(-0.005800 ,-0.000010 ,0.003200 ,0.000000),   array(-0.005900 ,0.000000 ,0.002600 ,0.000000),
        array(-0.005100 ,0.000000 ,0.002700 ,0.000000),    array(-0.003800 ,0.000000 ,0.001600 ,0.000000),  array( 0.002900 ,0.000000 ,-0.000100 ,0.000000),   array( 0.002900 ,0.000000 ,-0.001200 ,0.000000),
        array(-0.003100 ,0.000000 ,0.001300 ,0.000000),    array( 0.002600 ,0.000000 ,-0.000100 ,0.000000), array( 0.002100 ,0.000000 ,-0.001000 ,0.000000),   array( 0.001600 ,0.000000 ,-0.000800 ,0.000000),
        array(-0.001300 ,0.000000 ,0.000700 ,0.000000),    array(-0.001000 ,0.000000 ,0.000500 ,0.000000),  array(-0.000700 ,0.000000 ,0.000000 ,0.000000),    array( 0.000700 ,0.000000 ,-0.000300 ,0.000000),
        array(-0.000700 ,0.000000 ,0.000300 ,0.000000),    array(-0.000800 ,0.000000 ,0.000300 ,0.000000),  array( 0.000600 ,0.000000 ,0.000000 ,0.000000),    array( 0.000600 ,0.000000 ,-0.000300 ,0.000000),
        array(-0.000600 ,0.000000 ,0.000300 ,0.000000),    array(-0.000700 ,0.000000 ,0.000300 ,0.000000),  array( 0.000600 ,0.000000 ,-0.000300 ,0.000000),   array(-0.000500 ,0.000000 ,0.000300 ,0.000000),
        array( 0.000500 ,0.000000 ,0.000000 ,0.000000),    array(-0.000500 ,0.000000 ,0.000300 ,0.000000),  array(-0.000400 ,0.000000 ,0.000000 ,0.000000),    array( 0.000400 ,0.000000 ,0.000000 ,0.000000),
        array(-0.000400 ,0.000000 ,0.000000 ,0.000000),    array(-0.000300 ,0.000000 ,0.000000 ,0.000000),  array( 0.000300 ,0.000000 ,0.000000 ,0.000000),    array(-0.000300 ,0.000000 ,0.000100 ,0.000000),
        array(-0.000300 ,0.000000 ,0.000100 ,0.000000),    array(-0.000200 ,0.000000 ,0.000100 ,0.000000),  array(-0.000300 ,0.000000 ,0.000100 ,0.000000),    array(-0.000300 ,0.000000 ,0.000100 ,0.000000),
        array( 0.000200 ,0.000000 ,-0.000100 ,0.000000),   array(-0.000200 ,0.000000 ,0.000100 ,0.000000),  array( 0.000200 ,0.000000 ,-0.000100 ,0.000000),   array(-0.000200 ,0.000000 ,0.000100 ,0.000000),
        array( 0.000200 ,0.000000 ,0.000000 ,0.000000),    array( 0.000200 ,0.000000 ,-0.000100 ,0.000000), array( 0.000100 ,0.000000 ,-0.000100 ,0.000000),   array(-0.000100 ,0.000000 ,0.000000 ,0.000000),
        array( 0.000100 ,0.000000 ,-0.000100 ,0.000000),   array(-0.000200 ,0.000000 ,0.000100 ,0.000000),  array(-0.000100 ,0.000000 ,0.000000 ,0.000000),    array( 0.000100 ,0.000000 ,-0.000100 ,0.000000),
        array(-0.000100 ,0.000000 ,0.000100 ,0.000000),    array(-0.000100 ,0.000000 ,0.000100 ,0.000000),  array( 0.000100 ,0.000000 ,0.000000 ,0.000000),    array( 0.000100 ,0.000000 ,0.000000 ,0.000000),
        array( 0.000100 ,0.000000 ,-0.000100 ,0.000000),   array(-0.000100 ,0.000000 ,0.000000 ,0.000000),  array(-0.000100 ,0.000000 ,0.000000 ,0.000000),    array( 0.000100 ,0.000000 ,0.000000 ,0.000000),
        array( 0.000100 ,0.000000 ,0.000000 ,0.000000),    array(-0.000100 ,0.000000 ,0.000000 ,0.000000),  array( 0.000100 ,0.000000 ,0.000000 ,0.000000),    array( 0.000100 ,0.000000 ,0.000000 ,0.000000),
        array(-0.000100 ,0.000000 ,0.000000 ,0.000000),    array(-0.000100 ,0.000000 ,0.000000 ,0.000000),  array(-0.000100 ,0.000000 ,0.000000 ,0.000000),    array(-0.000100 ,0.000000 ,0.000000 ,0.000000),
        array(-0.000100 ,0.000000 ,0.000000 ,0.000000),    array(-0.000100 ,0.000000 ,0.000000 ,0.000000),  array(-0.000100 ,0.000000 ,0.000000 ,0.000000),    array( 0.000100 ,0.000000 ,0.000000 ,0.000000),    
        array(-0.000100 ,0.000000 ,0.000000 ,0.000000),    array( 0.000100 ,0.000000 ,0.000000 ,0.000000)
     );
     $secrad = M_PI/(180.0*3600.0);
     self::clcfundarg($tmjd,$fundarg);
     $dt = ($tmjd-51544.5e0)/36525.0e0;
     $deltapsi = 0.0;
     $deltaeps = 0.0;
     for (  $k = 0 ;  $k <= 105 ; $k++ ){
         $r = $nutarg[$k][0]*($fundarg[1]-$fundarg[4]);
         for (  $i = 1 ;  $i <= 4 ; $i++ )
           $r = $r + $nutarg[$k][$i] * $fundarg[$i+1];
         $deltapsi = $deltapsi+($coefnut[$k][0] + $dt*$coefnut[$k][1])*sin($r);
         $deltaeps = $deltaeps+($coefnut[$k][2] + $dt*$coefnut[$k][3])*cos($r);
       }
     $deltapsi = $secrad*($deltapsi); /* in radian */
     $deltaeps = $secrad*($deltaeps);
   }
  /* clcfundarg is a void
     to calculate fundamental arguments in radians
        on moment tmjd in modified julian date.
        fundarg[1..5] is massive of the values of arguments.
        global const - grarad - for change from degrees to radians. */
   static function clcfundarg ( $tmjd, $fundarg){
     $dt = ($tmjd-51544.5e0)/36525.0e0;
     $dt2 = $dt*$dt;
     $dt3 = $dt*$dt2;
     $fundarg[1] = 218.31643250e0+481267.8812772222e0*$dt-0.00161167e0*$dt2+0.00000528e0*$dt3;
     $fundarg[2] = 134.96298139e0+477198.8673980556e0*$dt+0.00869722e0*$dt2+0.00001778e0*$dt3;
     $fundarg[3] = 357.52772333e0+35999.05034e0*$dt-0.00016028e0*$dt2-0.00000333e0*$dt3;
     $fundarg[4] =  93.27191028e0+483202.0175380555e0*$dt-0.00368250e0*$dt2+0.00000306e0*$dt3;
     $fundarg[5] = 297.85036306e0+445267.11148e0*$dt-0.00191417e0*$dt2+0.00000528e0*$dt3;
     for (  $i = 1 ;  $i <= 5 ; $i++ )
       $fundarg[$i] = deg2rad($fundarg[$i]);
   }   
  static function clcprecangles ($epoch ,$tmjd ,&$dzita0 , &$teta ,&$zet ){
    $de = ($epoch-51544.5)/36525.0; /* in century */
    $de2 = $de*$de;
    $dt = ($tmjd-$epoch)/36525.0;    /* in century */
    $dt2 = $dt*$dt;
    $dt3 = $dt2*$dt;
    $r = (2306.2181+1.39656*$de-0.000139*$de2)*$dt; /* in second of arc */
    $dzita0 = $r+(0.30188-0.000344*$de)*$dt2+0.017998*$dt3;
    $teta = (2004.3109-0.85330*$de-0.000217*$de2)*$dt
          -(0.42665+0.000217*$de)*$dt2-0.041833*$dt3;
    $zet = $r+(1.09468+0.000066*$de)*$dt2+0.018203*$dt3;
  }
 static function descfromsphercoor($azt,$alt,$ro){
   $cb = cos($alt);
   $posdescart[1] = $ro*cos($azt)*$cb;
   $posdescart[2] = $ro*sin($azt)*$cb;
   $posdescart[3] = $ro*sin($alt);
   return $posdescart;
 }
 /* from descart position */
 static function sphericoorindegree ($poscart,$alpha,$delta,$range){
   $qq = ($poscart[1]*$poscart[1])+($poscart[2]*$poscart[2]);  
   $rr = $qq+($poscart[3]*$poscart[3]);
   $r = sqrt($rr);
   $q = sqrt($qq);
   $p = 1/$q;
   $f = 1/$r;
   $c = $p*$poscart[1];
   $s = $p*$poscart[2];
   $alpha = rad2deg(atan2(s,c));
   if($alpha < 0)$alpha += 360;
   $c = $q*$f;
   $s = $poscart[3]*$f;
   $delta = rad2deg(atan2($s,$c)); /* unforfun in degree */
   $range = $r;
 }
   
}