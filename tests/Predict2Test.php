<?php defined('SYSPATH') OR die('Kohana bootstrap needs to be included before tests run');

class Predict2Test extends Unittest_TestCase{

  function providerPredict2Test01(){
    $tle1 = array(
      'TEST SAT SGP 001',
      '1 88888U          80275.98708465  .00073094  13844-3  66816-4 0     9',
      '2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518   103',
      0
    );
    $tle2 = array(
      'TEST SAT SDP 001',
      '1 11801U          80230.29629788  .01431103  00000-0  14311-1 0     2',
      '2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848     2',
      64
    );
    return array(
      array(
          'tle'  => $tle1,'step' => 0.0,
          'x'    => 2328.97048951,'y'    => -5995.22076416,'z'    => 1719.97067261,
          'vx'   => 2.91207230,   'vy'   => -0.98341546,   'vz'   => -7.09081703
      ),array(
          'tle'  => $tle1,'step' => 360.0,
          'x'    => 2456.10705566,'y'    => -6071.93853760,'z'    => 1222.89727783,
          'vx'   => 2.67938992,   'vy'   => -0.44829041,   'vz'   => -7.22879231
      ),array(
          'tle'  => $tle1,'step' => 720.0,
          'x'    => 2567.56195068,'y'    => -6112.50384522,'z'    => 713.96397400,
          'vx'   => 2.44024599,   'vy'   => 0.09810869,    'vz'   => -7.31995916
      ),array(
          'tle'  => $tle1,'step' => 1080.0,
          'x'    => 2663.09078980,'y'    => -6115.48229980,'z'    => 196.39640427,
          'vx'   => 2.19611958,   'vy'   => 0.65241995,    'vz'   => -7.36282432
      ),array(
          'tle'  => $tle1,'step' => 1440.0,
          'x'    => 2742.55133057,'y'    => -6079.67144775,'z'    => -326.38095856,
          'vx'   => 1.94850229,   'vy'   => 1.21106251,    'vz'   => -7.35619372
      ),array(
          'tle'  => $tle2,'step' => 0.0,
          'x'    => 7473.37066650,'y'    => 428.95261765,'z'    => 5828.74786377,
          'vx'   => 5.1071513,    'vy'   => 6.44468284,  'vz'   => -0.18613096
      ),array(
          'tle'  => $tle2,'step' => 360.0,
          'x'    => -3305.22537232,'y'    => 32410.86328125,'z'    => -24697.17675781,
          'vx'   => -1.30113538,   'vy'   => -1.15131518,   'vz'   => -0.28333528
      ),array(
          'tle'  => $tle2,'step' => 720.0,
          'x'    => 14271.28759766,'y'    => 24110.46411133,'z'    => -4725.76837158,
          'vx'   => -0.32050445,   'vy'   => 2.67984074,    'vz'   => -2.08405289
      ),array(
          'tle'  => $tle2,'step' => 1080.0,
          'x'    => -9990.05883789,'y'    => 22717.35522461,'z'    => -23616.890662501,
          'vx'   => -1.01667246,   'vy'   => -2.29026759,   'vz'   => 0.72892364
      ),array(
          'tle'  => $tle2,'step' => 1440.0,
          'x'    => 9787.86975097, 'y'    => 33753.34667969,'z'    => -15030.81176758,
          'vx'   => -1.09425966,   'vy'   => 0.92358845,    'vz'   => -1.52230928
      )
    );
  }
  /**
   * @dataProvider providerPredict2Test01
   */
  function testPredict2Test01($tlef,$step,$x,$y,$z,$vx,$vy,$vz){
    //$file = file('tests/test-001.tle');
    $tle  = new Predict2_Model_TLE($tlef[0], $tlef[1], $tlef[2]);
    $sat  = new Predict2_Sat($tle);
    $sgpsdp  = new Predict2_SGPSDP();

    $sgpsdp->SGP4($sat, $step);
    Predict2_Math::Convert_Sat_State($sat->pos, $sat->vel);
    $this->assertSame($tlef[3],$sat->flags & Predict2_SGPSDP::DEEP_SPACE_EPHEM_FLAG);
    
    $this->assertEquals($x,$sat->pos->x,'x',1e-2);
    $this->assertEquals($y,$sat->pos->y,'y',1e-2);
    $this->assertEquals($z,$sat->pos->z,'z',1e-2);
    $this->assertEquals($vx,$sat->vel->x,'vx',1e-2);
    $this->assertEquals($vy,$sat->vel->y,'vy',1e-2);
    $this->assertEquals($vz,$sat->vel->z,'vz',1e-2);
  }
  
  public function providerFromCartToKepler()
  {
    return array(
      array(
        'pos'=>new Predict2_Model_Vector(-6576.0294136,3500.4872775,-2.2709637512e-06),
        'vel'=>new Predict2_Model_Vector(-0.005180327707,-0.056075328306,7.3184502353),
        'PosVar'=>new Predict2_Model_Vector(7457.6272607,0.0031609794801,89.593422823),
        'AngVar'=>new Predict2_Model_Vector(151.97317354,70.441171861,-70.100122333),
      ),
    );
  }
   /**
   * @dataProvider providerFromCartToKepler
   */
  function testFromCartToKepler($Pos_Test,$Vel_Test,$PosVar_Test,$AngVar_Test){
    $PosVar = new Predict2_Model_Vector();$AngVar = new Predict2_Model_Vector();
    Predict2_Conv::FromCartToKepler($Pos_Test,$Vel_Test,$PosVar,$AngVar);
    $this->assertEquals($AngVar->x,$AngVar_Test->x,'$AngVar->x',1e-3);//<=0;
    $this->assertEquals($AngVar->y,$AngVar_Test->y,'$AngVar->y',1e-3);//<=0;
    $this->assertEquals($AngVar->z,$AngVar_Test->z,'$AngVar->z',1e-3);//<=0;
    $this->assertEquals($PosVar->x,$PosVar_Test->x,'$PosVar->x',1e-3);//<=0;
    $this->assertEquals($PosVar->y,$PosVar_Test->y,'$PosVar->y',1e-3);//<=0;
    $this->assertEquals($PosVar->z,$PosVar_Test->z,'$PosVar->z',1e-3);//<=0;
  }
   /**
   * @dataProvider providerFromCartToKepler
   */
  function testFromKeplerToCart($Pos_Test,$Vel_Test,$PosVar_Test,$AngVar_Test){
    $Pos = new Predict2_Model_Vector();$Vel = new Predict2_Model_Vector();
    Predict2_Conv::FromKeplerToCart($PosVar_Test,$AngVar_Test,$Pos,$Vel);
    $this->assertEquals($Pos->x,$Pos_Test->x,'$Pos->x',1e-3);
    $this->assertEquals($Pos->y,$Pos_Test->y,'$Pos->y',1e-3);
    $this->assertEquals($Pos->z,$Pos_Test->z,'$Pos->z',1e-3);
    $this->assertEquals($Vel->x,$Vel_Test->x,'$Vel->x',1e-3);
    $this->assertEquals($Vel->y,$Vel_Test->y,'$Vel->y',1e-3);
    $this->assertEquals($Vel->z,$Vel_Test->z,'$Vel->z',1e-3);
  }
  function providerFromDatetoMJD(){
    return array(
    array(
      'civil time'=>'29-02-2036 02:45:00',
      'civil part'=>array(
        'h'=>2,'mi'=>45,'sec'=>0.0,
        'd'=>29,'m'=>2,'y'=>2036),
      'MJD in UTC' =>64752.11458333 
    ));
  }
  /**
  * @dataProvider providerFromDatetoMJD
  */
  function testFromDatetoMJD($ts,$tp,$jmd1){
    //forJob31 from civil date to MJD
    $t2 = strtotime($ts);
    $jmd2=Predict2_Time::unix2daynum($t2); //julian day
    $jmd2-=2400000.5;//MJD
    $this->assertEquals($jmd2,$jmd1,'',1e-8);
    $jmd = Predict2_Conv::FromDatetoMJD($tp['d'],$tp['m'],$tp['y'],$tp['h'],$tp['mi'],$tp['sec']);
    $this->assertEquals($jmd,$jmd1,'',1e-8);
    //$imd3 = math::Mdata($tp['y'],$tp['m'],$tp['d']);
  }
  /**
  * @dataProvider providerFromDatetoMJD
  */
  function testTransMJDtoDate($ts,$tp,$jmd){
    // 3.2 From MJD to civil date
    $t1 = strtotime($ts);
    $t2 = Predict2_Conv::TransMJDtoDate($jmd);
    $this->assertSame($t1,$t2);
    $t3=Predict2_Time::daynum2unix($jmd + 2400000.5);
    $this->assertSame($t1,$t3);
    $this->assertSame($ts,date('d-m-Y H:i:s',$t3));
  }  
  function providerForJob140(){
    return array(
      array(
        'tle' =>array('GORIZONT 25             ',
            '1 21922U 92017A   04167.63115611 -.00000281  00000-0  00000+0 0  9741', 
            '2 21922   7.6381  58.1852 0001479  67.2109 145.7651  1.00284792 44659',
        ),
        'year'=>4,
        'dayp'=>167.63115611,
        'dn'=>-0.00000281,
        'ai'=>7.6381,
        'au'=>58.1852,
        'ae'=>0.0001479,
        'ao'=>67.2109,
        'am'=>145.7651,
        'an'=>1.00284792,
        'te'=>53171.6311561100 - 0.125,  // UTC+3 to UTC
        'xx'=>array('x'=>690.108,'y'=>-42050.058,'z'=>-3051.419),
        'vv'=>array('x'=>3.054349,'y'=>0.074685,'z'=>-0.342793),
      )
    );
  }
/**
  * @dataProvider providerForJob140
  * @outputBuffering disabled
  */
  function testForJob140_descard($ftle,$year,$dayp,$dn,$ai,$au,$ae,$ao,$am,$an,$te,$xt,$vt){
    $tle = new Predict2_Model_TLE($ftle[0],$ftle[1],$ftle[2]);
    $xx  = new Predict2_Model_Vector();
    $vv  = new Predict2_Model_Vector();
    $tx  = Predict2_Conv::FromNoradToDesCart($tle,0,$xx,$vv);
    $this->assertSame($te,$tx);
    $this->assertEquals($xt['x'],$xx->x,'x',1e-3);
    $this->assertEquals($xt['y'],$xx->y,'y',1e-3);
    $this->assertEquals($xt['z'],$xx->z,'z',1e-3);
    $this->assertEquals($vt['x'],$vv->x,'vx',1e-3);
    $this->assertEquals($vt['y'],$vv->y,'vy',1e-3);
    $this->assertEquals($vt['z'],$vv->z,'vz',1e-3);
  }
  /**
  * @dataProvider providerForJob140
  */
  function testForJob140_norad($ftle,$year,$dayp,$dn,$ai,$au,$ae,$ao,$am,$an,$te,$xt,$vt){
    //Обратное преобразование
    $tle = new Predict2_Model_TLE($ftle[0],$ftle[1],$ftle[2]);
    $tle1= new Predict2_Model_TLE();
    $xx  = new Predict2_Model_Vector($xt['x'],$xt['y'],$xt['z']);
    $vv  = new Predict2_Model_Vector($vt['x'],$vt['y'],$vt['z']);
    Predict2_Conv::FromDesCartToNorad($te,$xx,$vv,$tle1); 
    $this->assertSame($tle->epoch_day,$tle1->epoch_day,'day');
    $this->assertSame($tle->epoch_fod,$tle1->epoch_fod,'fod');
    $this->assertSame($tle->epoch_year,$tle1->epoch_year,'year');
    $this->assertEquals($tle->eo,$tle1->eo,'eo',1e-3);
    $this->assertEquals($tle->xno,$tle1->xno,'xno',1e-3);
    $this->assertEquals($tle->xmo,$tle1->xmo,'xmo',1e-1);
    $this->assertEquals($tle->xnodeo,$tle1->xnodeo,'xnodeo',1e-3);
    $this->assertEquals($tle->omegao,$tle1->omegao,'omegao',1e-1);
    $this->assertEquals($tle->xndt2o,$tle1->xndt2o,'xndt2o',1e-3);
    $this->assertEquals($tle->xndd6o,$tle1->xndd6o,'xndd6o',1e-3);
  }
  
  function testSocketPredict(){
    $satellite = 'RS-15';
    $fp = stream_socket_client("udp://192.9.200.7:1211", $errno, $errstr);
    if (!$fp) {
        $err = "ERROR: $errno - $errstr<br />\n";
    } else {
        fwrite($fp, "GET_LIST\n");
        $output=explode("\n", fread($fp,1000));
        $this->assertContains($satellite,$output);

        fwrite($fp, "GET_TLE $satellite\n");
        $output=explode("\n", fread($fp,1000));
        $this->assertContains($satellite,$output);
        
        
        fwrite($fp, "GET_SAT $satellite\n");
        $output =explode("\n", fread($fp, 100));
        list($name, $lon, $lat, $az, $el, $aos_seconds, $foot) = $output;
        $this->assertSame($name,$satellite,'GET_SAT');
        
        $Commands[]="GET_DOPPLER";
        $Commands[]="RELOAD_TLE";
        $Commands[]="GET_SUN";
        $Commands[]="GET_MOON";
        $Commands[]="GET_MODE";
        $Commands[]="GET_QTH";
        
        fwrite($fp, "GET_TIME\$\n");
        $output =explode("\n", fread($fp, 100));
        $d = strtotime($output[2]);
        $this->assertSame(date('D M d H:i:s Y',$d),$output[2],'GET_TIME$');
        
        fwrite($fp, "GET_VERSION\n");
        $output =explode("\n", fread($fp, 100));
        $this->assertSame('2.2.3',$output[0],'GET_VERSION');
        
        fwrite($fp, "GET_TIME\n");
        $output =explode("\n", fread($fp, 100));
        
        fwrite($fp, "GET_SAT_POS $satellite\n");
        $output =explode("\n", fread($fp, 100));

        fwrite($fp, "PREDICT\n");
        $output =explode("\n", fread($fp, 100));
        
        fclose($fp);
    }
  }
}