<?php defined('SYSPATH') OR die('No direct access allowed.');
/**
 * @package     Predict2
 * @copyright   (c) 2012 Kopex
 */
class Predict2_Conv
{
  /**
  * from NORAD elements to DesCart position and velosity 
  * @param Predict2_Model_TLE $tle norad params
  * @param int $dift delta in day
  * @param Predict2_Model_Vector $xx   x, y, z  coordinate in km
  * @param Predict2_Model_Vector $vv  vx,vy,vz velosity in km/s
  * @return float moment vector in MJD
  */
  static function FromNoradToDesCart(Predict2_Model_TLE $tle,$dift,
                                    Predict2_Model_Vector $xx,Predict2_Model_Vector $vv){
   $te=self::FromDateToMJD(1,1,$tle->epoch_year,0,0,0.0)-1.0; //{ before 1 january of year }
   $te=$te+$tle->epoch_fod + $tle->epoch_day;        //{ date in modified julian day }

   $ao = $tle->omegao;$ai = $tle->xincl;$au = $tle->xnodeo;
   $ae = $tle->eo;$am = $tle->xmo;$anob = $tle->xno;$dnob = $tle->xndt2o;
   $c=array();
   self::ElemC($ai,$ae,$anob,$c);
   $ddl=M_PI*$dnob*$dift*$dift; // +1/2 dn/dt (t-t_0)^2 radian }
   $tc =$te + $dift; //in day
   $sdift=$dift*86400.0;//in second
   $udl=deg2rad($am)+$c[42]*$sdift+$ddl; // in radian }
   $udg=deg2rad($ao)+$c[43]*$sdift+M_PI_2;
   $udh=deg2rad($au)+$c[44]*$sdift-M_PI_2;
   $o=1.0e0;
   $w=2.00000e0;
   $ol=0.0e0;
   $fm=$c[1];
   $rz=$c[2];
   $cj2=$c[3];
   $sig=$c[4];
   $a1=$c[5];
   $a2=$c[6];
   $a3=$c[7];
   $ab=$c[8];
   $ae=$c[9];
   $asz=$c[10];
   $ac=$c[11];
   $ad=$c[12];
   $adz=$c[13];
   $ax=$c[14];
   $ag=$c[15];
   $at=$c[16];
   $ag1=$c[17];
   $ag2=$c[18];
   $an=$c[19];
   $aq=$c[20];
   $ac0=$c[21];
   $ac1=$c[22];
   $ac2=$c[23];
   $ac3=$c[24];
   $ac4=$c[25];
   $ak1=$c[26];
   $ak2=$c[27];
   $aq0=$c[28];
   $aq1=$c[29];
   $aq2=$c[30];
   $az=$c[31];
   $ap0=$c[32];
   $ap1=$c[33];
   $ap2=$c[34];
   $af0=$c[35];
   $af1=$c[36];
   $af2=$c[37];
   $af3=$c[38];
   $af4=$c[39];
   $alam=$c[40];
   $an0=$c[41];
   $oh= $c[45];
   $amu=$c[46];
   $akf2=$c[47];
   $akf4=$c[48];
   $akp2=$c[49];
   $akp4=$c[50];
   $rxx=1.0-$ax*$ax;
   $rx=Sqrt($rxx);
   $z2=$cj2*$cj2;
   $z0=$cj2*$sig;
   $sp=sin($udl);
   $se=$sp;
   $cp=cos($udl);
   $ce=$cp;
   $sf=$ol;
   $cf=$o;
   $ep=$udl;
   $ea=$ep;
   $pudl=$ep-$udl;
   $spp=$w*$sp*$cp;
   $cpp=$cp*$cp-$sp*$sp;
   $spr=$w*$spp*$cpp;
   $sff=$w*$sf*$cf;
   $j=0;
   $de=1.0e0;
   while ( ( $j < 15 ) and ( $de > 1.0e-15 ) ){
       $j=$j+1;
       $pudl=$ep-$udl;
       $spp=$w*$sp*$cp;
       $cpp=$cp*$cp-$sp*$sp;
       $spr=$w*$spp*$cpp;
       $sff=$w*$sf*$cf;
       $cff=$cf*$cf-$sf*$sf;
       $sfr=$w*$sff*$cff;
       $ee=$udl+$az*$se+$alam*$pudl-$ap1*$sp-$ap2*$spp
             -$af1*$sf+$af2*$sff+$af3*($sf*$cff+$cf*$sff)-$af4*$sfr;
       $se=sin($ee);
       $ce=cos($ee);
       $rb=$o/($o-$ax*$ce);
       $sp=$rx*$rb*$se;
       $cp=($ce-$ax)*$rb;
       $ep=$ee+atan(($sp*$ce-$cp*$se)/($cp*$ce+$sp*$se));
       $ef=$ep+$udg+$oh*$pudl+$akp2*$spp+$akp4*$spr-$akf2*$sff+$akf4*$sfr;
       $sf=sin($ef);
       $cf=cos($ef);
       $de=Abs($ea-$ee);
       $ea=$ee;
   }
   $bs=$ab*($o-$ae*$ce);
   $rt=$o/($o-$at*$cf);
   $bt=$rt*(-$asz*$cf+$ag);
   $r1=$bs*$bs+$z2;
   $r2=$o-$bt*$bt;
   $r3=Sqrt($r1*$r2);
   $r4=-$ac*$cf+$aq;
   $r6=$o/($r4*$r4+$sf*$sf);
   $r5=Sqrt($r6);
   $bw=$udh+1.570796326794896e0+$amu*$pudl+$aq1*$sf-$aq2*$sff;
   $bw=$bw+$ac1*$sp+$ac2*$spp+$ac3*($sp*$cpp+$cp*$spp)+$ac4*$spr;
   $sw=$r4*$r5;
   $cw=$sf*$r5;
   $su=sin($bw);
   $cu=cos($bw);
   $xx->x=$r3*($cw*$cu-$sw*$su);
   $xx->y=$r3*($cw*$su+$sw*$cu);
   $xx->z=$z0+$bs*$bt;
   $ri=$o/($bs*$bs+$z2*$bt*$bt);
   $r1=$o/$r1;
   $r2=$o/$r2;
   $q1=Sqrt($o-8*$ak2*$sp*$sp);
   $q2=($rx/($o+$ax*$cp));$q2*=$q2;
   $q3=$ab*$ae*$ag2;
   $vs=$q3*$ri*$q2*$sp*$q1;
   $q1=Sqrt($o-8*$ak1*$cf*$cf);
   $q2=($asz-$ag*$at)*$ag1;
   $vt=$q2*$ri*$rt*$rt*$sf*$q1;
   $vw=$a3*$r1*$r2;
   $r4=$bs*$vs*$r1-$bt*$vt*$r2;
   $vv->x=$xx->x*$r4-$xx->y*$vw;
   $vv->y=$xx->y*$r4+$xx->x*$vw;
   $vv->z=$bs*$vt+$vs*$bt;
   return $tc;
  }
  /**
  * from DesCart position velosity to NORAD elements
  * @param float $t      moment in modified julian day
  * @param predict2_model_vector $xx   x, y, z  coordinate in km
  * @param predict2_model_vector $vv  vx,vy,vz velosity in km/s
  * @param Predict2_Model_TLE $tle NORAD params
  */
  static function FromDesCartToNorad($te,Predict2_Model_Vector $xx,Predict2_Model_Vector $vv,
                                                        Predict2_Model_TLE $tle){
    $c=array();
    $fm= 3.98600448000e+5;
    $rz= 6.37813700000e+3;
    $cj2= 2.09728850857e+2;
    $sig=-3.55889244242e-02;
    $z2=$cj2*$cj2;
    $z0=$cj2*$sig;
    $x1=$xx->x;
    $x2=$xx->y;
    $x3=$xx->z-$z0;
    $v1=$vv->x;
    $v2=$vv->y;
    $v3=$vv->z;
    $zx=$x3*$x3;
    $r2=$x1*$x1+$x2*$x2+$zx;
    $u2=$v1*$v1+$v2*$v2+$v3*$v3;
    $r1=$x1*$v1+$x2*$v2+$x3*$v3;
    $f2=0.5e0*($r2-$z2+Sqrt(($r2-$z2)*($r2-$z2)+4.0e0*$z2*$zx));
    $bs=Sqrt($f2);
    $bt=$x3/$bs;
    $z4=$z2*$bt*$bt;
    $fi=2.0e0*$fm/($f2+$z4);
    $vt=$bs*$v3-$bt*$r1;
    $vs=$bs*$r1+$z2*$bt*$v3;
    $a1=$u2-$fi*($bs-$z0*$bt);
    $a3=$x1*$v2-$v1*$x2;
    $a2=$r1*$r1-$r2*$u2+$z2*$v3*$v3-$fi*($bs*$z4+$z0*$bt*$f2);
    $c[1]=$fm;  $c[2]=$rz;  $c[3]=$cj2;  $c[4]=$sig;
    $c[5]=$a1;  $c[6]=$a2;  $c[7]=$a3;
    self::Alei($c);
    self::Crona($c);
    $ab=$c[8];      $ae=$c[9];
    $asz=$c[10];    $ac=$c[11];
    $ad=$c[12];    $adz=$c[13];
    $ax=$c[14];     $ag=$c[15];
    $at=$c[16];    $ag1=$c[17];
    $ag2=$c[18];    $an=$c[19];
    $aq=$c[20];    $ac0=$c[21];
    $ac1=$c[22];   $ac2=$c[23];
    $ac3=$c[24];   $ac4=$c[25];
    $ak1=$c[26];   $ak2=$c[27];
    $aq0=$c[28];   $aq1=$c[29];
    $aq2=$c[30];    $az=$c[31];
    $ap0=$c[32];   $ap1=$c[33];
    $ap2=$c[34];   $af0=$c[35];
    $af1=$c[36];   $af2=$c[37];
    $af3=$c[38];   $af4=$c[39];
    $alam=$c[40];  $an0=$c[41];   $oh= $c[45];   $amu=$c[46];
    $akf2=$c[47];  $akf4=$c[48];  $akp2=$c[49];  $akp4=$c[50];
    $ce=(1.0e0-$bs/$ab)/$ae;
    $se=Sqrt(1.0e0-$ce*$ce);
    if($vs < 0.0)$se=-$se;
    $ea=ATan2($se,$ce);
    $cp=($ce-$ax)/(1.0e0-$ax*$ce);
    $sp=Sqrt(1.0e0-$cp*$cp);
    if($vs < 0.0)$sp=-$sp;
    $ep=ATan2($sp,$cp);
    $sr=($bt-$ag)/($asz-$bt*$at);
    $cr=Sqrt(1.0e0-$sr*$sr);
    if($vt < 0.0)$cr=-$cr;
    $sf=$cr;   $cf=-$sr;
    $ef=ATan2($sf,$cf);
    $spp=2.0*$sp*$cp;     $cpp=$cp*$cp-$sp*$sp;
    $spr=2.0*$spp*$cpp;   $sff=2.0*$sf*$cf;
    $cff=$cf*$cf-$sf*$sf;    $sfr=2.0*$sff*$cff;
    $am0=$ea-$az*$se-$alam*$ep+$ap1*$sp+$ap2*$spp
          +$af1*$sf-$af2*$sff-$af3*($sf*$cff+$cf*$sff)+$af4*$sfr;
    $amr=$am0/(1-$alam);
    $pudl=$ep-$amr;
    $ao0=$ef-$ep-$oh*$pudl-$akp2*$spp-$akp4*$spr+$akf2*$sff-$akf4*$sfr;
    $aor=$ao0-M_PI_2;
    $r1=-$ac*$cf+$aq;
    $bw=ATan2($x2,$x1)-ATan2($r1,$sf);
    $au0=$bw-M_PI_2-$amu*$pudl
          -$ac1*$sp-$ac2*$spp-$ac3*($spp*$cp+$cpp*$sp)
          -$ac4*$spr-$aq1*$sf+$aq2*$sff;
    $aur=$au0+M_PI_2;
    $tle->xmo    =rad2deg(Predict2_Math::FMod2p($amr)); 
    $tle->omegao =rad2deg(Predict2_Math::FMod2p($aor));
    $tle->xnodeo =rad2deg(Predict2_Math::FMod2p($aur));
    $tle->xincl  =rad2deg(ATan2($asz,$ac)); 
    $tle->xno    =($an0/(1-$alam))*rad2deg(8.64e4/360.0);
    $tle->xndt2o =0.0;
    $tle->eo     =$ae;

    $t = getdate(Predict2_Conv::TransMJDtoDate($te-0.125));
    $tc= Predict2_Conv::FromDateToMJD(1,1,$t['year'],0,0,0.0)-1.0; //{ before 1 january of year }
    $tle->epoch_year = (int)$t['year'];
    $tle->epoch_day = (int)floor($te-$tc);
    $tle->epoch_fod = $te - $tc - $tle->epoch_day;
  }  
  /**
  * to change from
  * Descart position and velosity to Kepler osculating elements
  *  
  * @param predict_vector $Pos   position  x,y,z in km
  * @param predict_vector $Vel   velosity  vx,vy,vz in km/s
  * @return predict_vector $PosVar 1,2,3 : a in km, e,  Inclination in degree
  * @return predict_vector $AngVar 1,2,3 : the ascending node in degree
  *              : the argument of perigei in degree
  *              : the mean anomaly in degree 
  */
  static function FromCartToKepler(Predict2_Model_Vector $Pos,Predict2_Model_Vector $Vel,
                                    Predict2_Model_Vector $PosVar,Predict2_Model_Vector $AngVar){
    $r=Sqrt(($Pos->x*$Pos->x)+($Pos->y*$Pos->y)+($Pos->z*$Pos->z)); // in km }
    $v2=($Vel->x*$Vel->x)+($Vel->y*$Vel->y)+($Vel->z*$Vel->z); // in (km/s)**2 }
    $q=predict2::ge/$r; // fm in (km**3/s**2)/km from UnL01010 for the Earth }
    $h=$v2-2*$q; // (km/s)^2 energy integral }
    $c1=$Pos->y*$Vel->z - $Pos->z*$Vel->y; // three momentum integrals }
    $c2=$Pos->z*$Vel->x - $Pos->x*$Vel->z;
    $c3=$Pos->x*$Vel->y - $Pos->y*$Vel->x;
    $c=Sqrt($c1*$c1+$c2*$c2+$c3*$c3);
    $l1=-$q*$Pos->x + $Vel->y*$c3 - $Vel->z*$c2; // three integrals }
    $l2=-$q*$Pos->y + $Vel->z*$c1 - $Vel->x*$c3;
    $l3=-$q*$Pos->z + $Vel->x*$c2 - $Vel->y*$c1;
    $l=Sqrt($l1*$l1+$l2*$l2+$l3*$l3);
    
    $a=-predict2::ge/$h; // in km semi-major axis }
    $e=$l/predict2::ge; // eccentricity }
    $p=$c*$c/predict2::ge;
    $ci=$c3/$c; $si=Sqrt(Abs(1-$ci*$ci));
   If($si < 1.0e-11){
       $si=1.0e-11; $rincl=0; // singularity inclination = 0 }
   }Else  $rincl=atan2($si,$ci); // inclination in radian }
   $sn=$c1/($c*$si); $cn=-$c2/($c*$si);
   $node=atan2($sn,$cn); // node in radian }
   $sp=$l3/($l*$si); $cp=($l1/$l)*$cn+($l2/$l)*$sn;
   $argp=atan2($sp,$cp); // argument of perigei in radian }
   $q=1/$r; $su=$q*$Pos->z/$si; $cu=$q*$Pos->x*$cn + $q*$Pos->y*$sn;
   $sv=$su*$cp-$cu*$sp; $cv=$cu*$cp+$su*$sp; // argument of latitude u = v + omega }
   $TrueAnom=atan2($sv,$cv); // true anomaly in radian }
   $q=1/(1+$e*$cv); $se=Sqrt(1-$e*$e)*$sv*$q; $ce=($cv+$e)*$q;
   $EccAnom=$TrueAnom+aTan(($se*$cv-$ce*$sv)/($ce*$cv+$se*$sv));
   $MeanAnom=$EccAnom-$e*$se; // M = E - e * Sin(E) in radian }
   $PosVar->x=$a; // in km }
   $PosVar->y=$e;
   $PosVar->z=rad2deg($rincl); // in degree }
   $AngVar->x=rad2deg(Predict2_Math::FMod2p($node)); // in degree }
   $AngVar->y=rad2deg(Predict2_Math::FMod2p($argp)); // in degree }
   $AngVar->z=rad2deg(Predict2_Math::FMod2p($MeanAnom)); // in degree }
  }
  
  /**
  * to change from
  * Kepler osculating elements to Descart position and velosity
  * 
  * @param mixed $PosVar 1,2,3 
  *             : a,bp - большая полусь in km, 
  *             : e - эксцентриситет,  
  *             : i - Inclination in degree
  * @param mixed $AngVar 1,2,3 
  *             : gd0 - the ascending node in degree
  *             : wp0 - the argument of perigei in degree
  *             : m0 - the mean anomaly un degree
  * @return predict_vector $Pos   position  x,y,z in km
  * @return predict_vector $Vel   velosity  vx,vy,vz in km/s
  */
  static function FromKeplerToCart(Predict2_Model_Vector $PosVar,Predict2_Model_Vector $AngVar,
                                    Predict2_Model_Vector $Pos,Predict2_Model_Vector $Vel){
    $ElemA=$PosVar->x;
    $ElemE=$PosVar->y;
    $CosInc=Cos(deg2rad($PosVar->z));
    $SinInc=Sin(deg2rad($PosVar->z));
    $MeanAnom=deg2rad($AngVar->z); 
    $EccAnom=self::KeplerEquation($ElemE,$MeanAnom);
    $SinE=Sin($EccAnom);
    $CosE=Cos($EccAnom);
    $SinV=Sqrt(1-$ElemE*$ElemE)*$SinE/(1-$ElemE*$CosE);
    $CosV=($CosE-$ElemE)/(1-$ElemE*$CosE);
    $TrueAnom=$EccAnom+atan(($SinV*$CosE-$CosV*$SinE)/($CosV*$CosE+$SinV*$SinE));
    $ArgLatitude=$TrueAnom+predict2::de2ra*$AngVar->y;// { in radian }
    $SinU=Sin($ArgLatitude);
    $CosU=Cos($ArgLatitude);
    $AscendingNode=predict2::de2ra*$AngVar->x;
    $SinNode=Sin($AscendingNode);
    $CosNode=Cos($AscendingNode);
    $r=$ElemA*(1-$ElemE*$CosE); //{ in km }
    $xr=+$CosU*$CosNode-$SinU*$SinNode*$CosInc;
    $yr=+$CosU*$SinNode+$SinU*$CosNode*$CosInc;
    $zr=+$SinU*$SinInc;
    $xv=-$SinU*$CosNode-$CosU*$SinNode*$CosInc;
    $yv=-$SinU*$SinNode+$CosU*$CosNode*$CosInc;
    $zv=+$CosU*$SinInc;
    $q=Sqrt(predict2::ge/($ElemA*(1-($ElemE*$ElemE))));
    $vr=$q*$ElemE*$SinV;
    $vn=$q*(1+$ElemE*$CosV);
    $Pos->x=$r*$xr;
    $Pos->y=$r*$yr;
    $Pos->z=$r*$zr;
    $Vel->x=$vr*$xr+$vn*$xv;
    $Vel->y=$vr*$yr+$vn*$yv;
    $Vel->z=$vr*$zr+$vn*$zv;
  }
  static Function KeplerEquation( $ElemE , $MeanAnom){
    $cura=$MeanAnom+$ElemE*Sin($MeanAnom);// { in radian }
    $difa=1.0;
    $iter=0;
    While  ( ( $difa > 1.0e-15 ) AND ( $iter < 241 ) ){
       $EccAnom=$cura-($cura-$ElemE*Sin($cura)-$MeanAnom)/(1-$ElemE*Cos($cura));
       $difa=Abs($cura-$EccAnom);
       $cura=$EccAnom;
       $iter++;
    }
    return $EccAnom; //{ eccentric anomaly in radian }
  }

  static function Crona ( array &$c){
     $fm=$c[1];
     $cj2=$c[3];
     $sig=$c[4];
     $a=$c[8];
     $e=$c[9];
     $s=$c[10];
     $ac=Sqrt(Abs(1.0e0-$s*$s));
     if($c[7] < 0.0)$ac=-$ac;
     $o=1.0e0;
     $o2=2.0;
     $o3=3.0;
     $o4=0.25e0;
     $o5=0.5e0;
     $e2=$e*$e;
     $re=$o-$e2;
     $e4=$e2*$e2;
     $s2=$s*$s;
     $rs=$o-$s2;
     $s4=$s2*$s2;
     $p=$cj2/($a*$re);
     $p2=$p*$p;
     $ps=$p*$sig;
     $z2=$sig*$sig;
     $ax=$e*($o+$p2*$re*($o-$o2*$s2+$p2*($o3
             -16.0*$s2+14.0*$s4-$o2*$e2*$rs*$rs)));
     $c[14]=$ax;
     $ag=-$ps*($o-$o2*$s2-$p2*($o3-12.0*$s2+10.0*$s4+$e2*($o-$o2*$s4)));
     $at=$ps*$s*($o-$p2*(5.0-6.0*$s2-$e2*($o-$o2*$s2)));
     $c[15]=$ag;
     $c[16]=$at;
     $ak1=0.125e0*$p2*$s2*($re+$z2-4.0*$p2*$re*$rs);
     $ak2=0.125e0*$p2*$e2*($s2-$p2*($o-10.0*$s2+11.0e0*$s4+$e2*$s4));
     $c[26]=$ak1;
     $c[27]=$ak2;
     $pm=Sqrt($fm*$a*$re);
     $p2=$o5*$p2;
     $p4=$o5*$p2*$p2;
     $ag1=$pm*($o+$p2*($rs*($o3+$e2)+$z2*(6.0-7.0*$s2))
               -$p4*$rs*(9.0+11.0e0*$s2+$e2*(6.0+34.0*$s2)+$e4*($o+$o3*$s2)));
     $ag2=$pm*($o-$p2*($o3-4.0*$s2-$e2)-$p4*(9.0-72.0*$s2
               +64.0*$s4+$e2*($o2-40.0*$s2+48.0*$s4)+$e4));
     $an=$ag1*($o+$o2*$ak2+9.0*$ak2*$ak2)/($ag2*($o+$o2*$ak1+9.0*$ak1*$ak1));
     $c[17]=$ag1;
     $c[18]=$ag2;
     $c[19]=$an;
     $ps=$ps*$s*$ac;
     $p2=$p*$p;
     $p4=$p2/16.0e0;
     $c[20]=$o2*$ps*($o-$p2*(4.0-5.0*$s2+$e2*$s2));
     $r1=$p2*$ac;
     $r2=$o2*$p4;
     $r3=4.0*$r2;
     $r4=$o4*$r1*$p4;
     $ac0=-$r1*($o+$o5*$e2+$p4*(24.0-56.0*$s2
                -4.0*$e2*($o+16.0*$s2)-$e4*($o2+$o3*$s2)));
     $c[21]=$ac0;
     $c[22]=-$o2*$r1*$e*($o+$r2*(4.0-28.0*$s2-$e2*(6.0+7.0*$s2)));
     $c[23]=-$o4*$r1*$e2*($o-$r3*(11.0e0+$e2*($o+$s2)));
     $c[24]=16.0*$r4*$e*$e2*($o2-$s2);
     $c[25]=$r4*$e4*($o2+$s2);
     $aq0=-$o5*$r1*($re+$o3*$z2-$o2*$p4*$re*(30.0-35.0*$s2+$e2*($o2+$o3*$s2)));
     $c[28]=$aq0;
     $c[29]=$p2*$ps*$re;
     $c[30]=6.0*$r4*$s2*$re*$re;
     $f2=$o-$r3*$re*(4.0-7.0*$s2)-$p2*$r2*$re
          *(20.0-136.0*$s2+113.0*$s4-$e2*(8.0+24.0*$s2-47.0*$s4));
     $f1=$p2*$re*($o2-$o3*$s2+$r3*($o2-22.0*$s2+19.0*$s4
                  -$e2*$s2*(10.0-13.0*$s2)));
     $f0=$o/($f1+$f2);
     $f0=$f0/Sqrt($o-$ax*$ax);
     $f0=$f0*$re*$re;
     $b0=$r3*(-$s2+$r3*(6.0-20.0*$s2+15.0*$s4-$e2*$s4));
     $b1=-$p2*$r3*$s2*($o2-$s2);
     $b2=$o3*$p2*$r2*$s4;
     $ap0=$f0*($b0+$b1+$b2*($o+$o5*$ax*$ax));
     $c[32]=$ap0;
     $c[33]=$ax*($b1+$o2*$b2)*$f0;
     $c[34]=$o4*$e2*$b2*$f0;
     $f0=$p2*$f0*$ag2/$ag1;
     $af0=$o5*$f0*$s2*($o+$o3*$ak1);
     $c[35]=$af0;
     $c[36]=-$f0*$s*($o2*$ag-1.5e0*$s*$at);
     $c[38]=-$f0*$s2*$at/6.0;
     $c[37]=-$f0*$s2*($o4+$ak1);
     $c[39]=0.125e0*$f0*$s2*$ak1;
     $c[40]=-$ap0-$an*$af0;
     $c[45]=$an-$o;
     $c[42]=$c[41]/($o-$c[40]);
     $c[43]=$c[45]*$c[42];
     $c[46]=$ac0+$an*$aq0;
     $c[44]=$c[42]*$c[46];
     $c[47]=$ak1*($o+4.0*$ak1);
     $c[48]=-0.75e0*$ak1*$ak1;
     $c[49]=-$an*$ak2*($o+4.0*$ak2);
     $c[50]=0.75e0*$an*$ak2*$ak2;
     $c[23]=$c[23]+$aq0*$c[49];
     $c[30]=$c[30]+$aq0*$c[47];
     $c[34]=$c[34]+$af0*$c[49];
     $c[37]=$c[37]+$af0*$c[47];
  }

  static function ElemC ( $ai, $ae, $anob, array  &$c){
   $fm= 3.98600448000e+5;
   $rz= 6.37813700000e+3;
   $cj2= 2.09728850857e+2;
   $sig=-3.55889244242e-02;
   $an0=2*M_PI*$anob/86400.0;
   $s=deg2rad($ai);
   $asz=sin($s);
   $ac=cos($s);
   $ad=$asz;
   $adz=-$asz;
   $a1=-exp(2.0e0/3.0e0*log($fm*$an0));
   $f0=-$fm/$a1;
   $ab=$f0;
   $fc=$cj2/$f0;
   $re=1.0-$ae*$ae;
   $f1=2.0-$re;
   $f2=$re*$re;
   for($i=0;$i<4;$i++){
     $pe=$cj2/($ab*$re);
     $p2=$pe*$pe;
     $p4=$p2*$p2;
     $ps=$pe*$sig;
     $d2=$ad*$ad;
     $rd=1.0-$d2;
     $d4=$d2*$d2;
     $rc=1.0e0-2.0e0*$ps*$ad-$p2*$d2*$re;
     $rr=1.0e0+2.0e0*$p2*$d2*$f1+$p4*$d4*$f2;
     $r1=$rc/$rr;
     $f3=1.0-$rd*$r1;
     $ab=$f0*(1.0e0-$p2*$re*$rd*$r1);
     $dz2=$adz*$adz;
     $rdz=1.0-$dz2;
     $dz4=$dz2*$dz2;
     $rcz=1.0e0-2.0e0*$ps*$adz-$p2*$dz2*$re;
     $rrz=1.0e0+2.0e0*$p2*$dz2*$f1+$p4*$dz4*$f2;
     $r2=$rcz/$rrz;
     $f4=1.0-$rdz*$r2;
     $f5=$f3+$d2*$adz*($fc*$pe*(2.0e0*$ad+$adz)+2.0e0*$ps);
     $f6= $f4+$dz2*$ad*($fc*$pe*(2.0e0*$adz+$ad)+2.0e0*$ps);
     $f7=Sqrt(Abs(($f5/$f6)*($f4/$f3)));
     $f9=$r1/$r2;
     $f7_1=(1.0e0+$f7)*$asz;
     $f8=2.0e0*$f7*$f7_1*$f9;
     $adz=(-1.0e0+$f9-$f9*$f7_1*$f7_1+(1.0e0-$f9*$f7*$f7)*$dz2)/$f8;
     $ad=$f7_1+$f7*$adz;
   }
   $r3=$fm*$ab*$re*$rd*(1.0+2.0*$p2*$f1+$p4*$f2)*$r1;
   $r4=$fm*$ab*$re*$rdz*(1.0+2.0*$p2*$f1+$p4*$f2)*$r2;
   $a3=Sqrt(Abs($r3));
   if($ac < 0.0)$a3=-$a3;
   $a2=-$fm*$ab*$re*(1.0+2.0*$p2*$f1*$rd*$r1+$p4*$f2*$rd*$r1);
   $c[8]=$ab;
   $c[9]=$ae;
   $c[10]=$asz;
   $c[11]=$ac;
   $c[5]=$a1;
   $c[6]=$a2;
   $c[7]=$a3;
   $c[12]=$ad;
   $c[13]=$adz;
   $c[1]=$fm;
   $c[2]=$rz;
   $c[3]=$cj2;
   $c[4]=$sig;
   $c[31]=-$c[8]*$c[9]*$a1/$fm;
   $c[41]=$an0;
   self::Crona($c);
  }

  static function Alei(array &$c){
   $o=1.0e0;
   $a1=$c[5];
   $a2=$c[6];
   $a3=$c[7];
   $fm=$c[1];
   $cj2=$c[3];
   $sig=$c[4];
   $ab=-$fm/$a1;
   $fa=$ab;
   $fe=-$a2/$fm;
   $re=$fe/$fa;
   $pa=$a3*$a3;
   $rd=-$pa/$a2;
   $rdz=$rd;
   $fd=$rd;
   $pa=$pa/$fm;
   for($i=0;$i<5;$i++){
       $pe=$cj2/($ab*$re);
       $p2=$pe*$pe;
       $ps=$pe*$sig;
       $d2=$o-$rd;
       $ad=Sqrt(Abs($d2));
       $p1=2.0*$p2*(2.0-$re);
       $p2=$p2*$re;
       $p4=$p2*$p2;
       $rr=$o+$p1*$d2+$p4*$d2*$d2;
       $rc=$o-2.0*$ps*$ad-$p2*$d2;
       $rr=$rc/$rr;
       $rf=($o+$p1+$p4)*$rr;
       $rt=$o+($rf-$rr)*$rd;
       $ab=$fa*($o-$p2*$rd*$rr);
       $re=$fe/($ab*$rt);
       $rd=$fd*$rt/$rf;
       $dz2=$o-$rdz;
       $adz=-Sqrt(Abs($dz2));
       $rrz=$o+$p1*$dz2+$p4*$dz2*$dz2;
       $rcz=$o-2.0*$ps*$adz-$p2*$dz2;
       $rdz=$pa*$rrz/($ab*$re*($o+$p1+$p4)*$rcz);
   }
   $e2=$o-$re;
   $c[8]=$ab;
   $c[9]=Sqrt(Abs($e2));
   $d2=$o-$rd;
   $c[12]=Sqrt(Abs($d2));
   $dz2=$o-$rdz;
   $c[13]=-Sqrt(Abs($dz2));
   $p2=$pe*$pe*(3.0-4.0*$d2+$e2);
   $c[10]=$c[12]+$ps*$rd*($o-3.0*$ps*$c[12]-$p2);
   $c[11]=Sqrt(Abs($o-$c[10]*$c[10]));
   if($a3 < 0.0)$c[11]=-$c[11];
   $rf=$o/$fm;
   $r=-$a1*$a1*$a1;
   $c[41]=Sqrt($r)*$rf;
   $c[31]=-$c[8]*$c[9]*$a1*$rf;
  }
  /**
  * UTC to MJD
  * 
  * @param mixed $day
  * @param mixed $month
  * @param mixed $year
  * @param mixed $hour
  * @param mixed $minute
  * @param mixed $second
  */
  static function FromDatetoMJD ($day,$month,$year,$hour,$minute,$second){
    if($year < 85)$year=2000+$year;
    if($year < 100)$year=1900+$year;
    $j1=$year-1900; $m1=$month-3;
    if($m1 < 0){$m1=$m1+12; $j1=$j1-1;}
    $q=(int)(0.5+30.6*$m1)+(int)($j1/4)+$day;
    $t=15078.0+365.0*$j1+$q;//in UTC
    $t=$t+($hour+($minute+$second/60.0)/60.0)/24.0;
    return $t;
  }
  /**
  * Convert MJD to Unix
  * 
  * @param mixed $Tmjd
  * @return timestamp
  */
  static function TransMJDtoDate ($Tmjd){
   $z=$Tmjd;//+0.125; // from UTC to Moscow decret time }
   $sp=86400.0*($z - (int)$z);
   $du=(int)$z; $r=$du-15078.0;  $n=(int)$r;
   $z=$r/1461.01; $nz=(int)$z; $n1=$n-1461*$nz;
   $r=$n1;  $z1=$r/365.25; $nz1=(int)$z1;
   $NYear=4*$nz+$nz1+1900;
   if($n1==1461){
      $NMonth=2;  $NDay=29;
   }else{
      $n2=$n1-365*$nz1; $r=$n2; $r=($r-0.5)/30.6; $m1=(int)$r;
      $NMonth=$m1+3;
      $NDay=$n2-(int)(30.6*$NMonth-91.3);
      if($NMonth > 12){
        $NMonth=$NMonth-12; $NYear=$NYear+1;
      }
   }
    $r=$sp/3600; $IntHour=(int)$r;
    $r=60*($r - (int)$r); $IntMin=(int)$r; $Second=60*($r - (int)$r);
    if($Second > 59.999){$Second=0; $IntMin++;}
    if($IntMin > 59){$IntMin=$IntMin-60; $IntHour++;}
    return gmmktime($IntHour,$IntMin,$Second,$NMonth,$NDay,$NYear);
  }
}