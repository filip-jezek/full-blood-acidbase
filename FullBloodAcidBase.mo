within ;
package FullBloodAcidBase
model FiggeFencl3
    "Base class for plasma acidbase after Figge and Fencl 3.0, missing inputs for SID, PCO2, total Pi and total albumin"
/*
Rem: Figge-Fencl Quantitative Physicochemical Model
Rem: of Human Acid-Base Physiology Version 3.0 (8 October, 2012; www.Figge-Fencl.org).
Rem:
Rem: Copyright 2003 - 2013 James J. Figge. Update published 28 April, 2013;
Rem: Update published 27 October, 2013.
Rem: 
Rem: The program may be downloaded free of charge for academic and educational use only.
Rem: This program is not intended for clinical use or for the care of human subjects in clinical trials.
Rem: This program applies to plasma-like solutions containing albumin. 
Rem: The program does not account for the contribution of plasma globulins, and has not been tested with clinical data; hence 
Rem: the program is not suitable for clinical use.


Implemented in Modelica by Filip Jezek, FEE CTU in Prague, 2016
*/
          // only constants here
  protected
  constant Real kw = 0.000000000000044;

  constant Real Kc1 = 0.0000000000244
      "Kc1 is derived from the parameters in the Henderson-Hasselbalch equation. pK = 6.1; a = 0.230 mM / kPa; 1 Torr = 0.13332236842105 kPa. The value of Kc1 is 2.44E-11 (Eq / L)^2 / Torr.";
  constant Real Kc2 = 0.00000000011
      "Kc2 is calculated from Harned and Scholes (1941) for 37 degrees C and ionic strength 0.15 M. The value of Kc2 is 5.5E-11 mol / L x 2 = 1.1E-10 Eq / L.";

  constant Real K1 = 0.0122
      "K1, K2, and K3 for the phosphoric acid - phosphate system are from Sendroy and Rem: Hastings (1927). pK1 = 1.915; pK2 = 6.66; pK3 = 11.78.";
  constant Real K2 = 0.000000219
      "K1, K2, and K3 for the phosphoric acid - phosphate system are from Sendroy and Hastings (1927). pK1 = 1.915; pK2 = 6.66; pK3 = 11.78.";
  constant Real K3 = 0.00000000000166
      "K1, K2, and K3 for the phosphoric acid - phosphate system are from Sendroy and Hastings (1927). pK1 = 1.915; pK2 = 6.66; pK3 = 11.78.";

//Rem: Enter desired values for SID, PCO2, [ Pi tot ], and [ Albumin ] in the next four lines.
  public
  input Real SID( unit = "meq/l") "Strong ion difference. Normal value 39";
  input Real pCO2(unit = "mmHg") "CO2 partial pressure. Normal value 40";
  input Real Pi( unit = "mmol/l") "Total phosphate. Normal value 1.15";
  input Real alb( unit="g/dl") "Albumin concentration. Normal value 4.4";

  Real H( unit="eq/l")= 10 ^ (-pH);
  Real pH(start = 10, unit="1");
  Real HO( unit="eq/l") = kw / H;
  Real HCO3( unit="eq/l") = Kc1 * pCO2 / H;
  Real CO3(  unit="eq/l")= Kc2 * HCO3 / H;

  protected
  Real FNX = K1 * H^2 + 2 * K1 * K2 * H + 3 * K1 * K2 * K3;
  Real FNY = H^3 + K1 * H^2 + K1 * K2 * H + K1 * K2 * K3;
  Real FNZ = FNX / FNY;
  public
  Real P( unit = "meq/l")= Pi * FNZ;
  Real Netcharge = SID + 1000 * (H - HO - HCO3 - CO3) - P;

  Real NB = 0.4 * (1 - (1 / (1 + 10 ^ (pH - 6.9))))
      "NB accounts for histidine pK shift due to the NB transition";

  // constant Real albuminResidues[:] = cat(1,{-1 /*cysteine */,-98/*glutamic acid*/,-18/*tyrosine*/,+24/*arginine */, /* lysine >>>*/ 2, 2, 2, 2, 1, 50} ,ones(16) /*histidine residues*/,/* amino terminus and carboxyl terminus*/{1, 1});
  protected
  Real albConversion = 1000 * 10 * alb / 66500;
  constant Real albuminResidues[:] = cat(1,{-1,              -98,                 -18,            +24,                              2, 2, 2, 2, 1, 50}, ones(16),                                                                {1, -1});
  // Real albuminPks[:] = {8.5 /* CYST*/,3.9 /* GLUT*/,11.7 /* TYR*/,12.5 /* ARG*/,/*LYS >>>*/5.8, 6.15, 7.51, 7.685, 7.86, 10.3,/*HIST>>>*/7.12 - NB, 7.22 - NB, 7.1 - NB, 7.49 - NB, 7.01 - NB, 7.31, 6.75, 6.36, 4.85, 5.76, 6.17, 6.73, 5.82, 5.1, 6.7, 6.2, 8/* amino terminus */,3.1 /*carboxyl terminus*/};
  Real albuminPks[:] = {8.5,          3.9,          11.7,         12.5,                    5.8, 6.15, 7.51, 7.685, 7.86, 10.3,           7.12 - NB, 7.22 - NB, 7.1 - NB, 7.49 - NB, 7.01 - NB, 7.31, 6.75, 6.36, 4.85, 5.76, 6.17, 6.73, 5.82, 5.1, 6.7, 6.2, 8,                    3.1};
  Real albChrg[n](each unit = "meq/l") "charge of albumin per unit";
  constant Integer n = size(albuminResidues, 1);
  public
  Real atch = sum(albChrg)*albConversion
      "albumin total charge. Normal value -12.2678";

equation
  Netcharge + atch = 0;

  for i in 1:n loop
    if albuminResidues[i] > 0 then
      // positive charge
      albChrg[i] = albuminResidues[i] * (1 / (1 + 10 ^ (pH - albuminPks[i])));
    else
      // negative charge
      albChrg[i] = albuminResidues[i] * (1 / (1 + 10 ^ (- pH + albuminPks[i])));
    end if;
  end for;

  /* 
  // debug
  alb2 = -1 * (1 / (1 + 10 ^ (-(pH - 8.5))))
  - 98 * (1 / (1 + 10 ^ (-(pH - 3.9))))
  - 18 * (1 / (1 + 10 ^ (-(pH - 11.7))))
  + 24 * (1 / (1 + 10 ^ (pH - 12.5)))
  + 2 * (1 / (1 + 10 ^ (pH - 5.8)))
  + 2 * (1 / (1 + 10 ^ (pH - 6.15)))
  + 2 * (1 / (1 + 10 ^ (pH - 7.51)))
  + 2 * (1 / (1 + 10 ^ (pH - 7.685)))
  + 1 * (1 / (1 + 10 ^ (pH - 7.86)))
  + 50 * (1 / (1 + 10 ^ (pH - 10.3)))
  + (1 / (1 + 10 ^ (pH - 7.12 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.22 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.1 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.49 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.01 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.31)))
  + (1 / (1 + 10 ^ (pH - 6.75)))
  + (1 / (1 + 10 ^ (pH - 6.36)))
  + (1 / (1 + 10 ^ (pH - 4.85)))
  + (1 / (1 + 10 ^ (pH - 5.76)))
  + (1 / (1 + 10 ^ (pH - 6.17)))
  + (1 / (1 + 10 ^ (pH - 6.73)))
  + (1 / (1 + 10 ^ (pH - 5.82)))
  + (1 / (1 + 10 ^ (pH - 5.1)))
  + (1 / (1 + 10 ^ (pH - 6.7)))
  + (1 / (1 + 10 ^ (pH - 6.2)))
  + (1 / (1 + 10 ^ (pH - 8)))
  - 1 * (1 / (1 + 10 ^ (-(pH - 3.1))));
*/
  annotation (                                 Documentation(info="<html>
<pre><font style=\"color: #006400; \">Rem:&nbsp;Figge-Fencl&nbsp;Quantitative&nbsp;Physicochemical&nbsp;Model</font>
<font style=\"color: #006400; \">Rem:&nbsp;of&nbsp;Human&nbsp;Acid-Base&nbsp;Physiology&nbsp;Version&nbsp;3.0&nbsp;(8&nbsp;October,&nbsp;2012;&nbsp;www.Figge-Fencl.org).</font>
<font style=\"color: #006400; \">Rem:</font>
<font style=\"color: #006400; \">Rem:&nbsp;Copyright&nbsp;2003&nbsp;-&nbsp;2013&nbsp;James&nbsp;J.&nbsp;Figge.&nbsp;Update&nbsp;published&nbsp;28&nbsp;April,&nbsp;2013;</font>
<font style=\"color: #006400; \">Rem:&nbsp;Update&nbsp;published&nbsp;27&nbsp;October,&nbsp;2013.</font>
<font style=\"color: #006400; \">Rem:&nbsp;</font>
<font style=\"color: #006400; \">Rem:&nbsp;The&nbsp;program&nbsp;may&nbsp;be&nbsp;downloaded&nbsp;free&nbsp;of&nbsp;charge&nbsp;for&nbsp;academic&nbsp;and&nbsp;educational&nbsp;use&nbsp;only.</font>
<font style=\"color: #006400; \">Rem:&nbsp;This&nbsp;program&nbsp;is&nbsp;not&nbsp;intended&nbsp;for&nbsp;clinical&nbsp;use&nbsp;or&nbsp;for&nbsp;the&nbsp;care&nbsp;of&nbsp;human&nbsp;subjects&nbsp;in&nbsp;clinical&nbsp;trials.</font>
<font style=\"color: #006400; \">Rem:&nbsp;This&nbsp;program&nbsp;applies&nbsp;to&nbsp;plasma-like&nbsp;solutions&nbsp;containing&nbsp;albumin.&nbsp;</font>
<font style=\"color: #006400; \">Rem:&nbsp;The&nbsp;program&nbsp;does&nbsp;not&nbsp;account&nbsp;for&nbsp;the&nbsp;contribution&nbsp;of&nbsp;plasma&nbsp;globulins,&nbsp;and&nbsp;has&nbsp;not&nbsp;been&nbsp;tested&nbsp;with&nbsp;clinical&nbsp;data;&nbsp;hence&nbsp;</font>
<font style=\"color: #006400; \">Rem:&nbsp;the&nbsp;program&nbsp;is&nbsp;not&nbsp;suitable&nbsp;for&nbsp;clinical&nbsp;use.</font>
<p><br><br><code><font style=\"color: #006400; \">I</font>mplemented in Modelica by Filip Jezek, FEE CTU in Prague, 2016</code></p>
</html>", revisions="<html>
<pre><font style=\"color: #006400; \">Filip Jezek, 2016</font></pre>
</html>"));
end FiggeFencl3;

model FiggeFencl3Detailed "Extension for detailed albumin balance"
  extends FiggeFencl3;
// protected
//   constant Real albuminResidues[:] = cat(1,{-1,              -98,                 -18,            +24,                              2, 2, 2, 2, 1, 50}, ones(16),                                                                {1, -1});
//   // Real albuminPks[:] = {8.5 /* CYST*/,3.9 /* GLUT*/,11.7 /* TYR*/,12.5 /* ARG*/,/*LYS >>>*/5.8, 6.15, 7.51, 7.685, 7.86, 10.3,/*HIST>>>*/7.12 - NB, 7.22 - NB, 7.1 - NB, 7.49 - NB, 7.01 - NB, 7.31, 6.75, 6.36, 4.85, 5.76, 6.17, 6.73, 5.82, 5.1, 6.7, 6.2, 8/* amino terminus */,3.1 /*carboxyl terminus*/};
//   Real albuminPks[:] = {8.5,          3.9,          11.7,         12.5,                    5.8, 6.15, 7.51, 7.685, 7.86, 10.3,           7.12 - NB, 7.22 - NB, 7.1 - NB, 7.49 - NB, 7.01 - NB, 7.31, 6.75, 6.36, 4.85, 5.76, 6.17, 6.73, 5.82, 5.1, 6.7, 6.2, 8,                    3.1};
//   Real albChrg[n](each unit = "meq/l") "charge of albumin per unit";
//   constant Integer n = size(albuminResidues, 1);
//   Real atch = sum(albChrg)*albConversion "albumin total charge. Normal value -12.2678";
  protected
  Real albPositivePart[n], albNegativePart[n];
  Real albTotalPlusPart[n], albTotalMinusPart[n];
  public
  Real HCO3mEqL = HCO3*1000;
  Real albHAPlus = albConversion*sum(albPositivePart) "A0 + H+ = HA+";
  Real albAMinus = -albConversion*sum(albNegativePart) "A- + H+ = HA0";
  Real albA0 = (ATotPlus - albHAPlus) "A0 + H+ = HA+";
  Real albHA0 = (ATotMinus - albAMinus) "HA0 = A- + H+";
  Real ATotPlus = albConversion*sum(albTotalPlusPart)
      "Part of albumin, which could be positive";
  Real ATotMinus = -albConversion*sum(albTotalMinusPart)
      "Part of albumin, which could be negative.";

  Real ach = sum(albChrg);
  Real atch0 = -12.2678 "to demonstrate the low buffer strength of albumin";
  Real atot1 = albHAPlus + albA0;
  Real atot2 = albAMinus + albHA0;
  Real test = atot1 + atot2;

  Real barGraphAlb1 = albHA0;
  Real barGraphAlb2 = albHA0 + albAMinus;
  Real barGraphAlb3 = albHA0 + albAMinus + albHAPlus;
  Real barGraphAlb4 = albHA0 + albAMinus + albHAPlus + albA0;
  Real barGraphHCO3Alb1 = P;
  Real barGraphHCO3Alb2 = P -atch;
  Real barGraphHCO3Alb3 = P -atch + 1000*HCO3;
  discrete Real dHCO3s, dHCO3e, dHCO3;
  discrete Real dAlbs, dAlbe, dAlb;
  Real diff = dHCO3 - dAlb;
  Real outofcobntrol = albHAPlus - albAMinus;
  parameter Real p = 1;
  parameter Real a = 1;
  parameter Real b = 1;
  parameter Real c = 1;
  Real AlbXMinus = -(0.148*pH - 0.818)*(alb*10) "Total albumin charge";
  Real PXminus = -(0.309*pH - 0.469)*Pi;
  Real Pminus = -P;
  Real totalDiff = AlbXMinus - atch + PXminus - Pminus;
/*
  Real dpH = der(pH);
  Real dAtch = der(atch)/dpH;
  Real dPi = der(P)/dpH;
  Real dbsHCO3 = der(HCO3*1000)/dpH;
  */
equation

  // POSSIBLE OPENMODELICA INCOMPATIBILITY
  // IN OPENMODELICA USE ELSEWHEN instead of end when; and when terminal.. lines
  when initial() then
    dHCO3s =  1000*HCO3;
    dAlbs = albHAPlus - albAMinus - p*P - 1000 * CO3*a - 1000 * HO*b + 1000*H*c;
  end when;
  when terminal() then
  // elsewhen terminal() then
    dHCO3e = 1000*HCO3;
    dHCO3 =  dHCO3e - dHCO3s;
    dAlbe = albHAPlus - albAMinus - p*P - 1000 * CO3*a - 1000 * HO*b + 1000*H*c;
    dAlb = dAlbe - dAlbs;
  end when;

  for i in 1:n loop
    if albuminResidues[i] > 0 then
      // positive charge
      albPositivePart[i] = albChrg[i];
      albNegativePart[i] = 0;
      albTotalPlusPart[i] = albuminResidues[i];
      albTotalMinusPart[i] = 0;
    else
      // negative charge
      albPositivePart[i] = 0;
      albNegativePart[i] = albChrg[i];
      albTotalPlusPart[i] = 0;
      albTotalMinusPart[i] = albuminResidues[i];
    end if;
  end for;

  /* 
  // debug
  alb2 = -1 * (1 / (1 + 10 ^ (-(pH - 8.5))))
  - 98 * (1 / (1 + 10 ^ (-(pH - 3.9))))
  - 18 * (1 / (1 + 10 ^ (-(pH - 11.7))))
  + 24 * (1 / (1 + 10 ^ (pH - 12.5)))
  + 2 * (1 / (1 + 10 ^ (pH - 5.8)))
  + 2 * (1 / (1 + 10 ^ (pH - 6.15)))
  + 2 * (1 / (1 + 10 ^ (pH - 7.51)))
  + 2 * (1 / (1 + 10 ^ (pH - 7.685)))
  + 1 * (1 / (1 + 10 ^ (pH - 7.86)))
  + 50 * (1 / (1 + 10 ^ (pH - 10.3)))
  + (1 / (1 + 10 ^ (pH - 7.12 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.22 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.1 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.49 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.01 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.31)))
  + (1 / (1 + 10 ^ (pH - 6.75)))
  + (1 / (1 + 10 ^ (pH - 6.36)))
  + (1 / (1 + 10 ^ (pH - 4.85)))
  + (1 / (1 + 10 ^ (pH - 5.76)))
  + (1 / (1 + 10 ^ (pH - 6.17)))
  + (1 / (1 + 10 ^ (pH - 6.73)))
  + (1 / (1 + 10 ^ (pH - 5.82)))
  + (1 / (1 + 10 ^ (pH - 5.1)))
  + (1 / (1 + 10 ^ (pH - 6.7)))
  + (1 / (1 + 10 ^ (pH - 6.2)))
  + (1 / (1 + 10 ^ (pH - 8)))
  - 1 * (1 / (1 + 10 ^ (-(pH - 3.1))));
*/
  annotation (                                 Documentation(info="<html>
<pre><font style=\"color: #006400; \">Rem:&nbsp;Figge-Fencl&nbsp;Quantitative&nbsp;Physicochemical&nbsp;Model</font>
<font style=\"color: #006400; \">Rem:&nbsp;of&nbsp;Human&nbsp;Acid-Base&nbsp;Physiology&nbsp;Version&nbsp;3.0&nbsp;(8&nbsp;October,&nbsp;2012;&nbsp;www.Figge-Fencl.org).</font>
<font style=\"color: #006400; \">Rem:</font>
<font style=\"color: #006400; \">Rem:&nbsp;Copyright&nbsp;2003&nbsp;-&nbsp;2013&nbsp;James&nbsp;J.&nbsp;Figge.&nbsp;Update&nbsp;published&nbsp;28&nbsp;April,&nbsp;2013;</font>
<font style=\"color: #006400; \">Rem:&nbsp;Update&nbsp;published&nbsp;27&nbsp;October,&nbsp;2013.</font>
<font style=\"color: #006400; \">Rem:&nbsp;</font>
<font style=\"color: #006400; \">Rem:&nbsp;The&nbsp;program&nbsp;may&nbsp;be&nbsp;downloaded&nbsp;free&nbsp;of&nbsp;charge&nbsp;for&nbsp;academic&nbsp;and&nbsp;educational&nbsp;use&nbsp;only.</font>
<font style=\"color: #006400; \">Rem:&nbsp;This&nbsp;program&nbsp;is&nbsp;not&nbsp;intended&nbsp;for&nbsp;clinical&nbsp;use&nbsp;or&nbsp;for&nbsp;the&nbsp;care&nbsp;of&nbsp;human&nbsp;subjects&nbsp;in&nbsp;clinical&nbsp;trials.</font>
<font style=\"color: #006400; \">Rem:&nbsp;This&nbsp;program&nbsp;applies&nbsp;to&nbsp;plasma-like&nbsp;solutions&nbsp;containing&nbsp;albumin.&nbsp;</font>
<font style=\"color: #006400; \">Rem:&nbsp;The&nbsp;program&nbsp;does&nbsp;not&nbsp;account&nbsp;for&nbsp;the&nbsp;contribution&nbsp;of&nbsp;plasma&nbsp;globulins,&nbsp;and&nbsp;has&nbsp;not&nbsp;been&nbsp;tested&nbsp;with&nbsp;clinical&nbsp;data;&nbsp;hence&nbsp;</font>
<font style=\"color: #006400; \">Rem:&nbsp;the&nbsp;program&nbsp;is&nbsp;not&nbsp;suitable&nbsp;for&nbsp;clinical&nbsp;use.</font>
<p><br><br><code><font style=\"color: #006400; \">I</font>mplemented in Modelica by Filip Jezek, FEE CTU in Prague, 2016</code></p>
</html>", revisions="<html>
<pre><font style=\"color: #006400; \">Filip Jezek, 2016</font></pre>
</html>"));
end FiggeFencl3Detailed;

  model FiggeFenclNSID "Calculation of normal SID using plasma figge fencl"
  //  extends FiggeFencl3Base;

    parameter Real pH0(unit = "mmHg") = 7.4;
    input Real pCO20(unit = "mmHg") = 40;
    input Real Pi0( unit = "mmol/l");//= 1.15;
    input Real alb0( unit="g/dl");//= 4.4;
    Real SID;// = figgeFencl3.SID;

    FiggeFencl3 figgeFencl3(pH=pH0, pCO2 = pCO20, Pi = Pi0, alb = alb0, SID = SID)
      annotation (Placement(transformation(extent={{-58,0},{-38,20}})));
  equation
  //  figgeFencl3Base.pH = pH0;
  //   pCO20 = figgeFencl3.pCO2;
  //   Pi0 = figgeFencl3.Pi;
  //   alb0 = figgeFencl3.alb;

  end FiggeFenclNSID;

  model FiggeFencl3pCO2 "Change of pH while changing pCO2"
   // extends FiggeFencl3Base;

    FiggeFencl3Detailed figgeFencl3Base(
      SID=SID,
      pCO2=pCO2,
      Pi=Pi,
      alb=alb) annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
    parameter Real SID = 39;
    Real pCO2 = time*40 + 20;
    parameter Real Pi = 1.15;
    parameter Real alb = 4.4;
    output Real pH = figgeFencl3Base.pH;
    annotation (experiment(Tolerance=1e-006), __Dymola_experimentSetupOutput);
  end FiggeFencl3pCO2;

model FiggeFencl3Extended
    "Base class for plasma acidbase after Figge and Fencl 3.0, missing inputs for SID, PCO2, total Pi and total albumin"
/*
Rem: Figge-Fencl Quantitative Physicochemical Model
Rem: of Human Acid-Base Physiology Version 3.0 (8 October, 2012; www.Figge-Fencl.org).
Rem:
Rem: Copyright 2003 - 2013 James J. Figge. Update published 28 April, 2013;
Rem: Update published 27 October, 2013.
Rem: 
Rem: The program may be downloaded free of charge for academic and educational use only.
Rem: This program is not intended for clinical use or for the care of human subjects in clinical trials.
Rem: This program applies to plasma-like solutions containing albumin. 
Rem: The program does not account for the contribution of plasma globulins, and has not been tested with clinical data; hence 
Rem: the program is not suitable for clinical use.


Implemented in Modelica by Filip Jezek, FEE CTU in Prague, 2016
*/
          // only constants here
  protected
  constant Real kw = 0.000000000000044;

  constant Real Kc1 = 0.0000000000244
      "Kc1 is derived from the parameters in the Henderson-Hasselbalch equation. pK = 6.1; a = 0.230 mM / kPa; 1 Torr = 0.13332236842105 kPa. The value of Kc1 is 2.44E-11 (Eq / L)^2 / Torr.";
  constant Real Kc2 = 0.00000000011
      "Kc2 is calculated from Harned and Scholes (1941) for 37 degrees C and ionic strength 0.15 M. The value of Kc2 is 5.5E-11 mol / L x 2 = 1.1E-10 Eq / L.";

  constant Real K1 = 0.0122
      "K1, K2, and K3 for the phosphoric acid - phosphate system are from Sendroy and Rem: Hastings (1927). pK1 = 1.915; pK2 = 6.66; pK3 = 11.78.";
  constant Real K2 = 0.000000219
      "K1, K2, and K3 for the phosphoric acid - phosphate system are from Sendroy and Hastings (1927). pK1 = 1.915; pK2 = 6.66; pK3 = 11.78.";
  constant Real K3 = 0.00000000000166
      "K1, K2, and K3 for the phosphoric acid - phosphate system are from Sendroy and Hastings (1927). pK1 = 1.915; pK2 = 6.66; pK3 = 11.78.";

//Rem: Enter desired values for SID, PCO2, [ Pi tot ], and [ Albumin ] in the next four lines.
  public
  input Real SID( unit = "meq/l") "Strong ion difference. Normal value 39";
  input Real pCO2(unit = "mmHg") "CO2 partial pressure. Normal value 40";
  input Real Pi( unit = "mmol/l") "Total phosphate. Normal value 1.15";
  input Real alb( unit="g/dl") "Albumin concentration. Normal value 4.4";

  Real H( unit="eq/l")= 10 ^ (-pH);
  Real pH(start = 10, unit="1");
  Real HO( unit="eq/l") = kw / H;
  Real HCO3( unit="eq/l") = Kc1 * pCO2 / H;
  Real CO3(  unit="eq/l")= Kc2 * HCO3 / H;

  protected
  Real FNX = K1 * H^2 + 2 * K1 * K2 * H + 3 * K1 * K2 * K3;
  Real FNY = H^3 + K1 * H^2 + K1 * K2 * H + K1 * K2 * K3;
  Real FNZ = FNX / FNY;
  public
  Real P( unit = "meq/l")= Pi * FNZ;
  Real Netcharge = SID + 1000 * (H - HO - HCO3 - CO3) - P -X;

  Real NB = 0.4 * (1 - (1 / (1 + 10 ^ (pH - 6.9))))
      "NB accounts for histidine pK shift due to the NB transition";

  // constant Real albuminResidues[:] = cat(1,{-1 /*cysteine */,-98/*glutamic acid*/,-18/*tyrosine*/,+24/*arginine */, /* lysine >>>*/ 2, 2, 2, 2, 1, 50} ,ones(16) /*histidine residues*/,/* amino terminus and carboxyl terminus*/{1, 1});
  protected
  Real albConversion = 1000 * 10 * alb / 66500;
  constant Real albuminResidues[:] = cat(1,{-1,              -98,                 -18,            +24,                              2, 2, 2, 2, 1, 50}, ones(16),                                                                {1, -1});
  // Real albuminPks[:] = {8.5 /* CYST*/,3.9 /* GLUT*/,11.7 /* TYR*/,12.5 /* ARG*/,/*LYS >>>*/5.8, 6.15, 7.51, 7.685, 7.86, 10.3,/*HIST>>>*/7.12 - NB, 7.22 - NB, 7.1 - NB, 7.49 - NB, 7.01 - NB, 7.31, 6.75, 6.36, 4.85, 5.76, 6.17, 6.73, 5.82, 5.1, 6.7, 6.2, 8/* amino terminus */,3.1 /*carboxyl terminus*/};
  Real albuminPks[:] = {8.5,          3.9,          11.7,         12.5,                    5.8, 6.15, 7.51, 7.685, 7.86, 10.3,           7.12 - NB, 7.22 - NB, 7.1 - NB, 7.49 - NB, 7.01 - NB, 7.31, 6.75, 6.36, 4.85, 5.76, 6.17, 6.73, 5.82, 5.1, 6.7, 6.2, 8,                    3.1};
  Real albChrg[n](each unit = "meq/l") "charge of albumin per unit";
  constant Integer n = size(albuminResidues, 1);
  public
  Real atch = sum(albChrg)*albConversion
      "albumin total charge. Normal value -12.2678";

    Real X = Xi*(Kx + Kxh*H)/H;
    parameter Real Kx = 0.00000000011;
    parameter Real Kxh = 0;
    parameter Real Xi = 1;
equation
  Netcharge + atch = 0;

  for i in 1:n loop
    if albuminResidues[i] > 0 then
      // positive charge
      albChrg[i] = albuminResidues[i] * (1 / (1 + 10 ^ (pH - albuminPks[i])));
    else
      // negative charge
      albChrg[i] = albuminResidues[i] * (1 / (1 + 10 ^ (- pH + albuminPks[i])));
    end if;
  end for;

  /* 
  // debug
  alb2 = -1 * (1 / (1 + 10 ^ (-(pH - 8.5))))
  - 98 * (1 / (1 + 10 ^ (-(pH - 3.9))))
  - 18 * (1 / (1 + 10 ^ (-(pH - 11.7))))
  + 24 * (1 / (1 + 10 ^ (pH - 12.5)))
  + 2 * (1 / (1 + 10 ^ (pH - 5.8)))
  + 2 * (1 / (1 + 10 ^ (pH - 6.15)))
  + 2 * (1 / (1 + 10 ^ (pH - 7.51)))
  + 2 * (1 / (1 + 10 ^ (pH - 7.685)))
  + 1 * (1 / (1 + 10 ^ (pH - 7.86)))
  + 50 * (1 / (1 + 10 ^ (pH - 10.3)))
  + (1 / (1 + 10 ^ (pH - 7.12 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.22 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.1 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.49 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.01 + NB)))
  + (1 / (1 + 10 ^ (pH - 7.31)))
  + (1 / (1 + 10 ^ (pH - 6.75)))
  + (1 / (1 + 10 ^ (pH - 6.36)))
  + (1 / (1 + 10 ^ (pH - 4.85)))
  + (1 / (1 + 10 ^ (pH - 5.76)))
  + (1 / (1 + 10 ^ (pH - 6.17)))
  + (1 / (1 + 10 ^ (pH - 6.73)))
  + (1 / (1 + 10 ^ (pH - 5.82)))
  + (1 / (1 + 10 ^ (pH - 5.1)))
  + (1 / (1 + 10 ^ (pH - 6.7)))
  + (1 / (1 + 10 ^ (pH - 6.2)))
  + (1 / (1 + 10 ^ (pH - 8)))
  - 1 * (1 / (1 + 10 ^ (-(pH - 3.1))));
*/
  annotation (                                 Documentation(info="<html>
<pre><font style=\"color: #006400; \">Rem:&nbsp;Figge-Fencl&nbsp;Quantitative&nbsp;Physicochemical&nbsp;Model</font>
<font style=\"color: #006400; \">Rem:&nbsp;of&nbsp;Human&nbsp;Acid-Base&nbsp;Physiology&nbsp;Version&nbsp;3.0&nbsp;(8&nbsp;October,&nbsp;2012;&nbsp;www.Figge-Fencl.org).</font>
<font style=\"color: #006400; \">Rem:</font>
<font style=\"color: #006400; \">Rem:&nbsp;Copyright&nbsp;2003&nbsp;-&nbsp;2013&nbsp;James&nbsp;J.&nbsp;Figge.&nbsp;Update&nbsp;published&nbsp;28&nbsp;April,&nbsp;2013;</font>
<font style=\"color: #006400; \">Rem:&nbsp;Update&nbsp;published&nbsp;27&nbsp;October,&nbsp;2013.</font>
<font style=\"color: #006400; \">Rem:&nbsp;</font>
<font style=\"color: #006400; \">Rem:&nbsp;The&nbsp;program&nbsp;may&nbsp;be&nbsp;downloaded&nbsp;free&nbsp;of&nbsp;charge&nbsp;for&nbsp;academic&nbsp;and&nbsp;educational&nbsp;use&nbsp;only.</font>
<font style=\"color: #006400; \">Rem:&nbsp;This&nbsp;program&nbsp;is&nbsp;not&nbsp;intended&nbsp;for&nbsp;clinical&nbsp;use&nbsp;or&nbsp;for&nbsp;the&nbsp;care&nbsp;of&nbsp;human&nbsp;subjects&nbsp;in&nbsp;clinical&nbsp;trials.</font>
<font style=\"color: #006400; \">Rem:&nbsp;This&nbsp;program&nbsp;applies&nbsp;to&nbsp;plasma-like&nbsp;solutions&nbsp;containing&nbsp;albumin.&nbsp;</font>
<font style=\"color: #006400; \">Rem:&nbsp;The&nbsp;program&nbsp;does&nbsp;not&nbsp;account&nbsp;for&nbsp;the&nbsp;contribution&nbsp;of&nbsp;plasma&nbsp;globulins,&nbsp;and&nbsp;has&nbsp;not&nbsp;been&nbsp;tested&nbsp;with&nbsp;clinical&nbsp;data;&nbsp;hence&nbsp;</font>
<font style=\"color: #006400; \">Rem:&nbsp;the&nbsp;program&nbsp;is&nbsp;not&nbsp;suitable&nbsp;for&nbsp;clinical&nbsp;use.</font>
<p><br><br><code><font style=\"color: #006400; \">I</font>mplemented in Modelica by Filip Jezek, FEE CTU in Prague, 2016</code></p>
</html>", revisions="<html>
<pre><font style=\"color: #006400; \">Filip Jezek, 2016</font></pre>
</html>"));
end FiggeFencl3Extended;

 package Thrash

   model FullBloodBase
      "Extension of Figge Fencl by erythrocytes from Siggard andersen"
   //
   // FULL BLOOD EXTENSION
   // Real Hb( unit = "g/dl") = 15 "Hemoglobin";
   // Real Hk( unit = "1")= Hb / 33.34 "Hematocrit";
   // Real sO2 "O2 saturation";
   //
   // Real BEox = BE  - 0.2*(1-sO2)*Hb
   //     "Oxygen Saturation correction by Siggaard-andersen";
   // Real BE( unit = "meq/l") "Base Excess, dependent to oxygen saturation";

     Real BEp( unit = "meq/l")= BE - mHCO3/(1-Hct);
     Real BEe( unit = "meq/l")= BEox + mHCO3/Hk
        "Base excess outside erythrocytes is lowered";
     Real NSID( unit = "meq/l")= normalPlasma.SID
        "Normal plasma SID for given Alb and Pi at standard pCO2 and pH";

     Real mHCO3( unit = "meq/l")
        "charge of HCO3 which moved to plasma from erythrocyte compartment";

   // Real pHz2;
   //
   // Real pHz1 = ((BEe - (-24.26 + hco3 -7.4 * (9.5 + 1.63 * HbM)) * (1.0 -0.0143 * HbM))/ ((1.0 -0.0143 * HbM) * (9.5 + 1.63 * HbM)));
   // Real hco3 = 0.0304*PCO2*10^(pHe-6.1);

   // protected
   //   constant Real HbM = 33.34;
   //   constant Real p[:] = {-1.159e-009, 1.328e-008, 2.228e-007, 1.479e-005,-0.0005606,0.04644,-2.431};
   //   constant Real ph[:] =  {8.229e-009,-8.913e-008,-1.82e-006,-0.0001034,0.003499,-0.3109,19.5915};
   //   Real k = p[1]*BEe^6 + p[2]*BEe^5 + p[3]*BEe^4 + p[4]*BEe^3 + p[5]*BEe^2 + p[6]*BEe + p[7];  // these equations are from Kofranek, 1980
   //   Real h = ph[1]*BEe^6 + ph[2]*BEe^5 + ph[3]*BEe^4 + ph[4]*BEe^3 + ph[5]*BEe^2 + ph[6]*BEe + ph[7];  // these equations are from Kofranek, 1980
   //   Real lpCO2 = log10(PCO2);//*Math.LOG10E;
   //   Real pHe( unit = "1")= (lpCO2-h)/k;

    protected
     .FullBloodAcidBase.FiggeFenclNSID normalPlasma(
       pH0=7.4,
       PCO20=40,
       alb0=alb,
       Pi0=Pi) annotation (Placement(transformation(extent={{-100,60},{-80,80}})));
    public
     FiggeFencl3 plasma(
       SID=SID,
       PCO2=pCO2,
       Pi=Pi,
       alb=alb) annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
     SAnomogram_formalization.SAoriginal fullErythrocyte
       annotation (Placement(transformation(extent={{-20,20},{0,40}})));

     parameter Real SID = 39;
     Real pCO2 = time*40 + 20;
     parameter Real Pi = 1.15;
     parameter Real alb = 4.4;
     output Real pH = fullErythrocyte.pH;

   equation
     plasma.pH = fullErythrocyte.pH;
     BEp = SID - NSID;
     annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
               -100},{100,100}}), graphics={
           Line(
             points={{-72,70},{-52,70},{-52,46}},
             color={0,0,255},
             smooth=Smooth.None,
             arrow={Arrow.None,Arrow.Filled}),
           Line(
             points={{-48,16},{-28,-38}},
             color={0,0,255},
             smooth=Smooth.None,
             arrow={Arrow.None,Arrow.Filled}),
           Line(
             points={{-8,16},{-24,-38}},
             color={0,0,255},
             smooth=Smooth.None,
             arrow={Arrow.None,Arrow.Filled})}));
   end FullBloodBase;

  model FullBloodAcidBase
    extends FullBloodBase;

  equation
    PCO2 = time*40 + 20;
    Pi = 1.15;
    alb = 4.4;
    BE = 0;
    Netcharge + atch = 0;
    sO2 = 1;
  end FullBloodAcidBase;

  model FullBloodpCO2
      "with rising pCO2, the amount of mHCO3 which moves to plasma rises, in excahnge to CL- (and thus mHCO3 charge decreases and SID RISES)"
    extends FullBloodBase;
  equation
    PCO2 = time*40 + 20;
    Pi = 1.15;
    alb = 4.4;
    BE = 0;
    Netcharge + atch = 0;
    sO2 = 1;
    annotation (experiment,
        __Dymola_experimentSetupOutput);
  end FullBloodpCO2;

  model FullBloodTitration
      "adding acid or base would change pH of blood. Compare with FF3 model wihtout erythrocytes"
    extends FullBloodBase;
  equation
    PCO2 = 40;
    Pi = 1.15;
    alb = 4.4;
    BE = time*60-30;
    Netcharge + atch = 0;
    sO2 = 1;
    annotation (experiment(
        StartTime=-20,
        StopTime=20,
        Tolerance=1e-006), __Dymola_experimentSetupOutput(derivatives=false));
  end FullBloodTitration;

  model FullBloodsO2 "low O2 saturation lowers BE"
    extends FullBloodBase;
  equation
    PCO2 = 40;
    Pi = 1.15;
    alb = 4.4;
    BE = 0;
    sO2 = time;
    Netcharge + atch = 0;
    annotation (experiment(StartTime=30, StopTime=50),
        __Dymola_experimentSetupOutput);
  end FullBloodsO2;

    model FiggeFencl3Graph
      extends FiggeFencl3;
    equation
      pH = time;
      HCO3 = Kc1 * PCO2 / H;
      annotation (experiment(StartTime=3, StopTime=9),
          __Dymola_experimentSetupOutput);
    end FiggeFencl3Graph;

    model FF3BicarbFlow
      extends FiggeFencl3;
      parameter Real bf = 0;
      Real bt;
    initial equation
      bt = 0;
    equation
      Netcharge + atch = 0;
      HCO3 = Kc1 * PCO2 / H + bt;
      der(bt) = bf;
    //  HCO3 = Kc1 * PCO2 / H;

    // NSID pro pCO 40 a ph 7.4
    // dSID = BEox = SID - NSID

    end FF3BicarbFlow;

    partial model KofrBase "Siggaard-andersen corrected to 37°C by Kofranek"

    constant Real Kw = 4.4E-14;
    constant Real Kc1 = 2.44E-11;
    constant Real Kc2 = 6.0E-11;
    constant Real K1 = 1.22E-2;
    constant Real K2 = 2.19E-7;
    constant Real K3 = 1.66E-12;
    Real glob = alb/1.5;
    Real H = 10^(-pH);
    Real Z = ( K1 * H*H          + 2* K1 * K2 * H + 3 * K1 * K2 * K3)  / ( H*H*H       + K1 * H*H + K1 * K2 * H + K1 * K2 * K3);
    Real NB = 0.4 * ( 1 - 1 / ( 1 + 10^( pH - 6.9)));
    Real DSIDglob = -0.75*glob*(pH-7.4);

    Real pH;
    Real SID;
    Real alb;
    Real Pi;
    Real pCO2;
    equation

      // upgraded to figge-fencl3.0 albumin charge
    SID = -(+DSIDglob + 1000 * ( H - Kw / H - Kc1 * pCO2 / H - Kc1 * Kc2 * pCO2 /(H^2))  - Pi * Z
     + ( -1 / ( 1 + 10^( - ( pH - 8.5)))
     - 98 / ( 1 + 10^( - ( pH - 3.9)))
     - 18 / ( 1 + 10^( - ( pH - 11.7)))
     + 24 / ( 1 + 10^( + ( pH - 12.5)))
     + 2 / ( 1 + 10^( pH - 5.80))
     + 2 / ( 1 + 10^( pH - 6.15))
     + 2 / ( 1 + 10^( pH - 7.51))
     + 2 / ( 1 + 10^( pH - 7.685))
     + 1 / ( 1 + 10^( pH - 7.86))
     + 50 / ( 1 + 10^( pH - 10.3))
     + 1 / ( 1 + 10^( pH - 7.12 + NB))
     + 1 / ( 1 + 10^( pH - 7.22 + NB))
     + 1 / ( 1 + 10^( pH - 7.1 + NB))
     + 1 / ( 1 + 10^( pH - 7.49 + NB))
     + 1 / ( 1 + 10^( pH - 7.01 + NB))
     + 1 / ( 1 + 10^( pH - 7.31))
     + 1 / ( 1 + 10^( pH - 6.75))
     + 1 / ( 1 + 10^( pH - 6.36))
     + 1 / ( 1 + 10^( pH - 4.85))
     + 1 / ( 1 + 10^( pH - 5.76))
     + 1 / ( 1 + 10^( pH - 6.17))
     + 1 / ( 1 + 10^( pH - 6.73))
     + 1 / ( 1 + 10^( pH - 5.82))
     + 1 / ( 1 + 10^( pH - 5.1))
     + 1 / ( 1 + 10^( pH - 6.7))
     + 1 / ( 1 + 10^( pH - 5.50))
     + 1 / ( 1 + 10^( pH - 8.0))
     - 1 / ( 1 + 10^( - ( pH - 3.1))))    * 1000 * 10 *  alb  / 66500);
                                       // from 4.0
                                  // from 6.0
                                  // from 7.6
                                   //from 7.8
                                  // from 8
                                       // from 7.19
                                       // from 7.29
                                      // from 7.17
                                       // from 7.56
                                       // from 7.08
                                       // from 7.38
                                       // from 6.82
                                       // from 6.43
                                       // from 4.92
                                       // from 5.83
                                       // from 6.24
                                       // from 6.80
                                       // from 5.89
                                      // from 5.20
                                      // from 6.80
                                       // from 5.50
                                       // from 8.0

    // // ORIGINAL COMPUTATION IN kOFRANEK FLASH SIMULATOR
    //   SID = -(+DSIDglob + 1000 * ( H - Kw / H - Kc1 * pCO2 / H - Kc1 * Kc2 * pCO2 /(H^2))  - Pi * Z
    //  + ( -1 / ( 1 + 10^( - ( pH - 8.5)))
    //   - 98 / ( 1 + 10^( - ( pH - 4.0)))
    //   - 18 / ( 1 + 10^( - ( pH - 11.7)))
    //   + 24 / ( 1 + 10^( + ( pH - 12.5)))
    //   + 2 / ( 1 + 10^( pH - 5.80))
    //   + 2 / ( 1 + 10^( pH - 6.00))
    //   + 1 / ( 1 + 10^( pH - 7.60))
    //   + 2 / ( 1 + 10^( pH - 7.80))
    //   + 2 / ( 1 + 10^( pH - 8.00))
    //   + 50 / ( 1 + 10^( pH - 10.3))
    //   + 1 / ( 1 + 10^( pH - 7.19 + NB))
    //   + 1 / ( 1 + 10^( pH - 7.29 + NB))
    //   + 1 / ( 1 + 10^( pH - 7.17 + NB))
    //   + 1 / ( 1 + 10^( pH - 7.56 + NB))
    //   + 1 / ( 1 + 10^( pH - 7.08 + NB))
    //   + 1 / ( 1 + 10^( pH - 7.38))
    //   + 1 / ( 1 + 10^( pH - 6.82))
    //   + 1 / ( 1 + 10^( pH - 6.43))
    //   + 1 / ( 1 + 10^( pH - 4.92))
    //   + 1 / ( 1 + 10^( pH - 5.83))
    //   + 1 / ( 1 + 10^( pH - 6.24))
    //   + 1 / ( 1 + 10^( pH - 6.80))
    //   + 1 / ( 1 + 10^( pH - 5.89))
    //   + 1 / ( 1 + 10^( pH - 5.20))
    //   + 1 / ( 1 + 10^( pH - 6.80))
    //   + 1 / ( 1 + 10^( pH - 5.50))
    //   + 1 / ( 1 + 10^( pH - 8.0))
    //   - 1 / ( 1 + 10^( - ( pH - 3.1))))    * 1000 * 10 *  alb  / 66500);
    end KofrBase;

    model KofrpCO2
      extends Thrash.KofrBase;

    equation
      SID = 38.5;
      alb = 4.4;
      Pi = 1.15;
      pCO2 = time;

      annotation (experiment(
          StartTime=30,
          StopTime=50,
          Tolerance=1e-006), __Dymola_experimentSetupOutput(derivatives=false));
    end KofrpCO2;

    model zander
    Real BE, Hb, hco3, pH, pH2, sO2, pCO2;
    equation
    BE = (1-0.0143*Hb)*(
        (hco3 - 24.26) +
        (1.63*Hb + 9.5)*(pH-7.4))
        -0.2*Hb*(1-sO2);
    Hb = time*15;
    sO2 = 1;
    //BE = time*40 - 20;
    BE = time*20 -10;
    //hco3 = 24.26;
    hco3 = time*24.6 + 24.26/2;
    //hco3 = 0.0304*pCO2*10^(pH-6.1);

    pCO2 = time*20 + 30;

    pH2 = ((BE - (-24.26 + hco3 -7.4 * (9.5 + 1.63 * Hb)) * (1.0 -0.0143 * Hb))/ ((1.0 -0.0143 * Hb) * (9.5 + 1.63 * Hb)));

    //pH2 = (BE - (-24.26 + hco3 -7.4 * (9.5 + 1.63 * Hb) -0.2 * Hb * (1.0 - sO2)))/(9.5 + 1.63 * Hb);

    end zander;

    partial model KofrBase1980 "model from publication, but not working"
    protected
     Real a[:] = {996.35 - 10.35*Hb,
     35.16875 + 0.25875*Hb,
     -82.41 + 2.01*Hb,
     -5.27625 - (5.025e-2)*Hb,
     121 - Hb,
     2.625 + 0.025*Hb,
     -2.556 - 0.0944*Hb,
     13.87634 + 0.186653*Hb + 0.00534936*Hb^2,
     0.548 - 0.0274*Hb,
     0.274 - 0.0137*Hb};

      constant Real p[:] = {-1.159e-009, 1.328e-008, 2.228e-007, 1.479e-005,-0.0005606,0.04644,-2.431};
      constant Real ph[:] =  {8.229e-009,-8.913e-008,-1.82e-006,-0.0001034,0.003499,-0.3109,19.5915};

      Real k = p[1]*BEe^6 + p[2]*BEe^5 + p[3]*BEe^4 + p[4]*BEe^3 + p[5]*BEe^2 + p[6]*BEe + p[7];  // these equations are from Kofranek (?)
      Real h = ph[1]*BEe^6 + ph[2]*BEe^5 + ph[3]*BEe^4 + ph[4]*BEe^3 + ph[5]*BEe^2 + ph[6]*BEe + ph[7];  // these equations are from Kofranek (?)
      Real lpCO2 = log10(pCO2);//*Math.LOG10E;
    public
     Real pHe( unit = "1")= (lpCO2-h)/k;
     Real BEe;

     Real BE = BEox + 0.2*(1-sO2)*Hb;
     Real Y = a[7] + sqrt(a[8] + a[9]*BE);
     Real pH_k1980 = (a[1]*a[10] + a[2]*Y + (a[3]*a[10] + a[4]*Y)*log10(pCO2))/
                 (a[5]*a[10]  +a[6]*Y);
     parameter Real sO2 = 1;
     Real Hct = Hb/33.35;
     Real Hb;// = 15;

     Real BEox;
     Real pCO2;

     Real pHz1, pHz2;

     Real hco3 = 0.0304*pCO2*10^(pHz1-6.1);

      // Zander 1995
     Real BE_z1 = (1-0.0143*Hb)*(
        (hco3 - 24.26) +
        (1.63*Hb + 9.5)*(pHz1-7.4))
        -0.2*Hb*(1-sO2);

      // zander lang 2002 (cited as zander 1995)
     Real BE_z2 = (1 - 0.0143*Hb)*(
                  (0.0304*pCO2*(10^(pHz2-6.1)) - 24.26)
                  + (9.5 + 1.63*Hb)*(pHz2 - 7.4))
                  - 0.2*Hb*(1-sO2);

    equation
      BE_z1 = BE_z2;
      BE_z1 = BE;
      BEe = BE;
      //BEox = time*40-20;//0;

      annotation (experiment,
          __Dymola_experimentSetupOutput);
    end KofrBase1980;

    model KofrComparisonTitration
      extends KofrBase1980;
    equation

        Hct = 15/33.34;//21;//time*33;
    //   pCO2 = time*40 + 20;
       pCO2 = 40;
    //    BEox = 0;
       BEox = time*60 - 30;

    end KofrComparisonTitration;

    model KofrComparisonpCO2
      extends KofrBase1980;
    equation

       Hct = 1;//21;//time*33;
       pCO2 = time*40 + 20;
       // pCO2 = 40;
       BEox = 0;
       // BEox = time*60 - 30;

    end KofrComparisonpCO2;

    partial model pCO2range_base

      parameter Real SID = 39;
      Real pCO2 = time*40 + 20;
      parameter Real Pi = 1.15;
      parameter Real alb = 4.4;
      output Real pH;
    equation

    end pCO2range_base;

    model FiggeFenclTitration "Change of pH while changing SID (or BE)"
      extends FiggeFencl3;
      Real NSID( unit = "meq3333/l")= normalPlasma.figgeFencl3Base.SID
        "Normal plasma SID for given Alb and Pi at standard pCO2 and pH";

          Real BEp = SID - NSID;
    .FullBloodAcidBase.FiggeFenclNSID normalPlasma(
        pH0=7.4,
        PCO20=40,
        alb0=alb,
        Pi0=Pi);
    equation
      PCO2 = 40;
      Pi = 1.15;
      alb = 4.4;

      Netcharge + atch = 0;
      BEp = time*60-30;
    //  sO2 = 1;

    end FiggeFenclTitration;

    package Interfaces
      type ion = enumeration(
          Na "  1",
          K "   2",
          Ca "  3",
          Cl "  4",
          HCO3 "5",
          Alb " 6",
          Pi "  7",
          XA "  8");
        connector Ionogram
          Real c[ion]( each unit = "mmol/l");
          flow Real f[ion]( each unit = "mol/s");
        end Ionogram;

      model bloodInterface

        Ionogram ig
          annotation (Placement(transformation(extent={{-16,-4},{4,16}})));
          Real sO2 = 1; // TODO
          Real Hb = 15;

          Real SID = ig.c[ion.Na] + ig.c[ion.K] + ig.c[ion.Ca] - ig.c[ion.Cl] - ig.c[ion.XA];
          Real NSID = 39;//normalPlasma.SID;
          Real BE = SID - NSID;
          Real BEox = BE - 0.2*(1-sO2)*Hb;
          //FiggeFencl3.FiggeFenclNSID normalPlasma(pH0 = 7.4, PCO20 = 40, alb0 = ig.c[ion.Alb], Pi0 = ig.c[ion.Pi]);
      equation
      //  SID =

        // electroneutrality must hold
        //ig.c[ion.Na] + ig.c[ion.K] + ig.c[ion.Ca] - ig.c[ion.Cl] - ig.c[ion.HCO3] - ig.c[ion.Alb] - ig.c[ion.Pi] - ig.c[ion.XA] = 0;

        for i in 1:size(ion, 1) loop
          // TODO per DISTRIBUTION VOLUME!!!
          ig.f[i] = der(ig.c[i]);
        end for;
      end bloodInterface;
    end Interfaces;

    model Full_Blood
      "Exact formalization of original Siggaard-andersen nomogram. Same as FiggeFencl3.SAnomogram_formalization.SAoriginal model, but with inputs instead of parameters."

    /*protected 
  constant Real pco2BBCoef[:] = {2.1125e-009,  -640.9926e-009,     72.7649e-006,     -3.2862e-003,    -38.1749e-003,      8.2352e+000,    -97.0551e+000};
  constant Real pco2BECoef[:] = {8.3975e-009,   -513.9503e-009,      3.8105e-006,    231.6888e-006,    -46.5581e-003,    353.7105e-003,     39.9871e+000};
  constant Real pHBBCoef[:] = {40.8936e-012,    -13.0063e-009,      1.6780e-006,   -111.7919e-006,      4.0776e-003,    -67.8274e-003,      7.2888e+000};
  constant Real pHBECoef[:] = {131.3315e-009,      2.5027e-006,    175.6144e-006,     11.9273e-003,      7.4001e+000};
*/
    protected
      constant Real pco2BBCoef[:] = {2.1125e-009,  -640.9926e-009,     72.7649e-006,     -3.2862e-003,    -38.1749e-003,      8.2352e+000,    -97.0551e+000};
      constant Real pco2BECoef[:] = {8.3975e-009,   -513.9503e-009,      3.8105e-006,    231.6888e-006,    -46.5581e-003,    353.7105e-003,     39.9871e+000};
      constant Real pHBBCoef[:] = {40.8936e-012,    -13.0063e-009,      1.6780e-006,   -111.7919e-006,      4.0776e-003,    -67.8274e-003,      7.2888e+000};
      constant Real pHBECoef[:] = {131.3315e-009,      2.5027e-006,    175.6144e-006,     11.9273e-003,      7.4001e+000};

    public
      Real pCO2BB( start = 96) = pco2BBCoef[1]*BB^6 + pco2BBCoef[2]*BB^5 + pco2BBCoef[3]*BB^4 + pco2BBCoef[4]*BB^3 + pco2BBCoef[5]*BB^2  + pco2BBCoef[6]*BB + pco2BBCoef[7];
      Real pCO2BE( start = 40) = pco2BECoef[1]*BE^6 + pco2BECoef[2]*BE^5 + pco2BECoef[3]*BE^4 + pco2BECoef[4]*BE^3 + pco2BECoef[5]*BE^2  + pco2BECoef[6]*BE + pco2BECoef[7];
      Real pHBB( start = 7) = pHBBCoef[1]*BB^6 + pHBBCoef[2]*BB^5 + pHBBCoef[3]*BB^4 + pHBBCoef[4]*BB^3 + pHBBCoef[5]*BB^2  + pHBBCoef[6]*BB + pHBBCoef[7];
      Real pHBE( start = 7) = pHBECoef[1]*BE^4 + pHBECoef[2]*BE^3 + pHBECoef[3]*BE^2 + pHBECoef[4]*BE + pHBECoef[5];

      Real BB = BE+0.42*Hb+41.7;
      Real BE = BEox+0.2*(1-sO2)*Hb;

      Real pH( start = 7.4) = (log10(pCO2) - log10(pCO2BB))*(pHBB - pHBE)/(log10(pCO2BB) - log10(pCO2BE))+ pHBB;
    /*
public 
  Real pCO2BB( start = 96) = pco2BBCoef[1]*BB^6 + pco2BBCoef[2]*BB^5 + pco2BBCoef[3]*BB^4 + pco2BBCoef[4]*BB^3 + pco2BBCoef[5]*BB^2  + pco2BBCoef[6]*BB + pco2BBCoef[7];
  Real pCO2BE( start = 40) = pco2BECoef[1]*BE^6 + pco2BECoef[2]*BE^5 + pco2BECoef[3]*BE^4 + pco2BECoef[4]*BE^3 + pco2BECoef[5]*BE^2  + pco2BECoef[6]*BE + pco2BECoef[7];
  Real pHBB( start = 7) = pHBBCoef[1]*BB^6 + pHBBCoef[2]*BB^5 + pHBBCoef[3]*BB^4 + pHBBCoef[4]*BB^3 + pHBBCoef[5]*BB^2  + pHBBCoef[6]*BB + pHBBCoef[7];
  Real pHBE( start = 7) = pHBECoef[1]*BE^4 + pHBECoef[2]*BE^3 + pHBECoef[3]*BE^2 + pHBECoef[4]*BE + pHBECoef[5];

  Real BB = BE+0.42*cHb+41.7;
  Real BE = BEox+0.2*(1-sO2)*cHb;
  Real pH( start = 7.4) = (log10(pCO2) - log10(pCO2BB))*(pHBB - pHBE)/(log10(pCO2BB) - log10(pCO2BE))+ pHBB;
*/
      Real Hb(unit="g/dl") = Hct*33.34;
      input Real BEox;

      input Real Hct;// = 20;
      parameter Real sO2 = 1;

      input Real pCO2;
    equation
    //    pCO2 = time*40 +20;
    //    BEox = 0;

    //  BB = time*100;
    //  BE = time*44-22;
    end Full_Blood;

   model SA_comparison_BE
     parameter Real Hct = 15/33.34;
     Real BEox = time*60 - 30;
     parameter Real sO2 = 1;
     parameter Real pCO2 = 40;
     SAnomogram_formalization.Zander1995 zander1995(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
     SAnomogram_formalization.SAoriginal sAoriginal(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{66,38},{86,58}})));
     SAnomogram_formalization.SAVanSlyke sAVanSlyke(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-14,40},{6,60}})));
   SAnomogram_formalization.Kofr2009 kofr2009(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-22,40},{-2,60}})));
     SAnomogram_formalization.SAVanSlyke77 sAVanSlyke77(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
   end SA_comparison_BE;

   model SA_comparison_plasma_pCO2
      "pH dependent on varying pCO2 by diffrent approximations. We take our implementation of SA nomogram as reference. Compare the object's pH"
     parameter Real Hct = 0;
     parameter Real BEox = 0;
     parameter Real sO2 = 1;
     Real pCO2 = time*40 + 20;
     SAnomogram_formalization.Zander1995 zander1995(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
     SAnomogram_formalization.Kofr2009 kofr2009(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-20,40},{0,60}})));
     SAnomogram_formalization.SAoriginal sAoriginal(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{34,40},{54,60}})));
     SAnomogram_formalization.SAVanSlyke sAVanSlyke(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{80,40},{100,60}})));
     SAnomogram_formalization.SAVanSlyke77 sAVanSlyke77(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
     FiggeFencl3 figgeFencl3(
       SID=SID,
       pCO2=pCO2,
       Pi=Pi,
       alb=alb) annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
     parameter Real SID = 39 + BEox;
     parameter Real Pi = 1.15;
     parameter Real alb = 4.4;
     annotation (experiment(Tolerance=0.001), __Dymola_experimentSetupOutput,
       __Dymola_Commands(file="def.mos" "def"));
   end SA_comparison_plasma_pCO2;

   model full_blood_test_pco2_at_BE
      import FullBloodAcidBase;
     constant Real fullHb = 33.34;
     parameter Real Hb = 15;
     parameter Real Hct = Hb/fullHb;

     FullBloodAcidBase.Full_Blood.Full_blood_combo fb(
        BE=BE,
        pCO2=pCO2,
        Pi=Pi,
        alb=alb)
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));

     FiggeFencl3 figgeFencl3Base(
       SID=SID,
       pCO2=pCO2,
       Pi=Pi,
       alb=alb) annotation (Placement(transformation(extent={{-60,20},{-40,40}})));

     FullBloodAcidBase.SAnomogram_formalization.SAVanSlyke sAVanSlyke(  pCO2 = pCO2, BEox = BE, Hct = Hct, Alb = alb)
       annotation (Placement(transformation(extent={{-18,40},{2,60}})));

     Real BE = 0;
     Real pCO2;
     //parameter Real pCO2 = 40;// time*40 + 20;
     Real Pi = 1.15;
     Real alb = 4.4;//time*4.4*2;
     output Real pH = fb.pH;

     Real SID;
    protected
     FiggeFenclNSID normalPlasma(
       pH0=7.4,
       pCO20=40,
       alb0=alb,
       Pi0=Pi)
       annotation (Placement(transformation(extent={{-40,60},{-20,80}})));

   // public
   //   discrete Real bepop( start = -20);
     discrete Real sTime( start = 0);
   //   parameter Real diff = 10;
   equation
     SID - normalPlasma.SID = BE;

     pCO2 = (time - sTime)*40 + 20;
   //   BE = bepop;
   //
   //   when sample(1, 1) then
   //      bepop = pre(bepop) + diff;
   //      sTime = time;
   //   end when;

   end full_blood_test_pco2_at_BE;

   model test_combo_chng_Hb
      "Test Figge-fencl plasma and SA full hemoatocrite blood during variable pCO2 against "

     FiggeFencl3Detailed                         plasma(
       SID=SID,
       pCO2=pCO2,
       Pi=Pi,
       alb=alb) "Plasma compartment (Hct = 0)"
       annotation (Placement(transformation(extent={{-80,8},{-60,28}})));
     SAnomogram_formalization.SAoriginal fullErythrocyte(
       Hct=1,
       BEox=BEe,
       pCO2=pCO2) "Full erythrocyte blood with Hct = 1"
       annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

     constant Real fullHb = 33.34;
     Real Hb = (time - sTime)*29 + 1;
     Real Hct = Hb/fullHb;
     parameter Real k = 1
        "pisvajcova konstanta - kolik mHCO3 zvostane v erytrocytu? Nebo jina nepresnost? Ideal 0.77";

     Real BEp( unit = "meq/l") = BE - mHCO3/(1-Hct);
     Real BEe( unit = "meq/l")= BE + k*mHCO3/Hct;
     Real BE = 0;//time*40 - 20;
     Real mHCO3;

     Real SID;
     Real pCO2 = dpic;
     Real Pi = 1.15;
     parameter Real alb = 4.4;
     output Real pH = fullErythrocyte.pH;

    FiggeFencl3Detailed                         plasma_only(
       SID=normalPlasma.SID + BE,
       pCO2=pCO2,
       Pi=Pi,
       alb=alb) "Plasma model acc to FiggeFencl for comparison only" annotation (Placement(transformation(extent={{60,-20},{80,0}})));

     SAnomogram_formalization.SAoriginal normalBloodSA(
       Hct=Hct,
       BEox=BE,
       pCO2=pCO2) "SA original for comparison only"   annotation (Placement(transformation(extent={{60,20},{80,40}})));

    protected
     FiggeFenclNSID normalPlasma(
       pH0=7.4,
       pCO20=40,
       alb0=alb,
       Pi0=Pi)
       annotation (Placement(transformation(extent={{-80,40},{-60,60}})));

   discrete Real dpic( start = 20);
   discrete Real sTime( start = 0);
   equation
     plasma.pH = fullErythrocyte.pH;
     BEp = SID - normalPlasma.SID;

     when sample(1,1) then
       sTime = time;
       dpic = pre(dpic) + 10;
     end when;
     annotation (experiment(
         StopTime=1,
         __Dymola_NumberOfIntervals=500,
         Tolerance=1e-003), Diagram(coordinateSystem(preserveAspectRatio=false,
             extent={{-100,-100},{100,100}}), graphics={
           Line(
             points={{-84,52},{-90,52},{-98,52},{-94,20},{-86,18}},
             color={28,108,200},
             arrow={Arrow.None,Arrow.Filled}),
           Line(
             points={{-54,16},{-30,2},{-26,-6}},
             color={28,108,200},
             arrow={Arrow.None,Arrow.Filled}),
           Line(
             points={{-54,-12},{-36,-14},{-26,-10}},
             color={28,108,200},
             arrow={Arrow.None,Arrow.Filled}),
           Text(
             extent={{-24,-14},{-4,0}},
             lineColor={28,108,200},
             fillColor={0,0,255},
             fillPattern=FillPattern.Solid,
             textString="pH"),
           Line(
             points={{10,84},{10,-74},{10,-80},{10,-82}},
             color={28,108,200},
             arrow={Arrow.None,Arrow.Filled})}));
   end test_combo_chng_Hb;

   model full_blood_test_BE_at_pco2
      import FullBloodAcidBase;
     constant Real fullHb = 33.34;
     parameter Real Hb = 15;
     parameter Real Hct = Hb/fullHb;

     FullBloodAcidBase.Full_Blood.Full_blood_combo fb(
        BE=BE,
        pCO2=pCO2,
        Pi=Pi,
        alb=alb)
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));

     FiggeFencl3 figgeFencl3Base(
       SID=SID,
       pCO2=pCO2,
       Pi=Pi,
       alb=alb) annotation (Placement(transformation(extent={{-60,20},{-40,40}})));

     FullBloodAcidBase.SAnomogram_formalization.SAVanSlyke sAVanSlyke(  pCO2 = pCO2, BEox = BE, Hct = Hct, Alb = alb)
       annotation (Placement(transformation(extent={{-18,40},{2,60}})));

     Real BE; //= 0;// time*40 - 20;
     Real pCO2;
     //parameter Real pCO2 = 40;// time*40 + 20;
     Real Pi = 1.15;
     Real alb = 4.4;//time*4.4*2;
     output Real pH = fb.pH;

     Real SID;
    protected
     FiggeFenclNSID normalPlasma(
       pH0=7.4,
       pCO20=40,
       alb0=alb,
       Pi0=Pi)
       annotation (Placement(transformation(extent={{-40,60},{-20,80}})));

    public
     discrete Real bepop( start = 20);
     discrete Real sTime( start = 0);
     parameter Real diff = 10;
   equation
     SID - normalPlasma.SID = BE;

     BE = (time - sTime)*10 - 20;
     pCO2 = bepop;

     when sample(1, 1) then
        bepop = pre(bepop) + diff;
        sTime = time;
     end when;

   end full_blood_test_BE_at_pco2;

   model full_blood_test_VS_alb_at_BE
      import FullBloodAcidBase;
     constant Real fullHb = 33.34;
     parameter Real Hb = 15;
     parameter Real Hct = Hb/fullHb;

     FullBloodAcidBase.Full_Blood.Full_blood_combo fb(
        BE=BE,
        pCO2=pCO2,
        Pi=Pi,
        alb=alb)
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));

     FiggeFencl3 figgeFencl3Base(
       SID=SID,
       pCO2=pCO2,
       Pi=Pi,
       alb=alb) annotation (Placement(transformation(extent={{-60,20},{-40,40}})));

     FullBloodAcidBase.SAnomogram_formalization.SAVanSlyke sAVanSlyke(  pCO2 = pCO2, BEox = BE, Hct = Hct, Alb = alb)
       annotation (Placement(transformation(extent={{-18,40},{2,60}})));

         FullBloodAcidBase.SAnomogram_formalization.SAVanSlyke77 sAVanSlyke77( pCO2 = pCO2, BEox = BE, Hct = Hct, sO2 = 1)
       annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
     Real BE; //= 0;// time*40 - 20;
     Real pCO2;
     //parameter Real pCO2 = 40;// time*40 + 20;
     Real Pi = 1.15;
     Real alb; //= 4.4;//time*4.4*2;
     output Real pH = fb.pH;

     Real SID;
    protected
     FiggeFenclNSID normalPlasma(
       pH0=7.4,
       pCO20=40,
       alb0=alb,
       Pi0=Pi)
       annotation (Placement(transformation(extent={{-40,60},{-20,80}})));

    public
     discrete Real bepop( start = - 20);
     discrete Real sTime( start = 0);
     parameter Real diff = 10;
   equation
     SID - normalPlasma.SID = BE;

     pCO2 = 40;
     // pCO2 = (time - sTime)*40 + 20;

     BE = bepop;
     alb = (time-sTime)*4.4*2;
     when sample(1, 1) then
        bepop = pre(bepop) + diff;
        sTime = time;
     end when;

   end full_blood_test_VS_alb_at_BE;

   model full_blood_test_VS_Pi_at_BE
      import FullBloodAcidBase;
     constant Real fullHb = 33.34;
     parameter Real Hb = 15;
     parameter Real Hct = Hb/fullHb;

     FullBloodAcidBase.Full_Blood.Full_blood_combo fb(
        BE=BE,
        pCO2=pCO2,
        Pi=Pi,
        alb=alb)
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));

     FiggeFencl3 figgeFencl3Base(
       SID=SID,
       pCO2=pCO2,
       Pi=Pi,
       alb=alb) annotation (Placement(transformation(extent={{-60,20},{-40,40}})));

     FullBloodAcidBase.SAnomogram_formalization.SAVanSlyke sAVanSlyke(  pCO2 = pCO2, BEox = BE, Hct = Hct, Alb = alb)
       annotation (Placement(transformation(extent={{-18,40},{2,60}})));

         FullBloodAcidBase.SAnomogram_formalization.SAVanSlyke77 sAVanSlyke77( pCO2 = pCO2, BEox = BE, Hct = Hct, sO2 = 1)
       annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
     Real BE; //= 0;// time*40 - 20;
     Real pCO2;
     //parameter Real pCO2 = 40;// time*40 + 20;
     Real Pi;  //= 1.15;
     Real alb = 4.4;//time*4.4*2;
     output Real pH = fb.pH;

     Real SID;
    protected
     FiggeFenclNSID normalPlasma(
       pH0=7.4,
       pCO20=40,
       alb0=alb,
       Pi0=Pi)
       annotation (Placement(transformation(extent={{-40,60},{-20,80}})));

    public
     discrete Real bepop( start = - 20);
     discrete Real sTime( start = 0);
     parameter Real diff = 10;
   equation
     SID - normalPlasma.SID = BE;

     pCO2 = 40;
     // pCO2 = (time - sTime)*40 + 20;

     BE = bepop;
     Pi = (time-sTime)*1.15*2;
     when sample(1, 1) then
        bepop = pre(bepop) + diff;
        sTime = time;
     end when;

     annotation (experiment(StopTime=5, Tolerance=1e-005),
         __Dymola_experimentSetupOutput(doublePrecision=true));
   end full_blood_test_VS_Pi_at_BE;

   model test_combo_chng_Alb
      "Test Figge-fencl plasma and SA full hemoatocrite blood during variable pCO2 against "

     FiggeFencl3                         plasma(
       SID=SID,
       pCO2=pCO2,
       Pi=Pi,
       alb=alb) "Plasma compartment (Hct = 0)"
       annotation (Placement(transformation(extent={{-80,8},{-60,28}})));
     SAnomogram_formalization.SAoriginal fullErythrocyte(
       Hct=1,
       BEox=BEe,
       pCO2=pCO2) "Full erythrocyte blood with Hct = 1"
       annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

     constant Real fullHb = 33.34;
     Real Hb = 15;
     Real Hct = Hb/fullHb;
     parameter Real k = 1
        "pisvajcova konstanta - kolik mHCO3 zvostane v erytrocytu? Nebo jina nepresnost? Ideal 0.77";

     Real BEp( unit = "meq/l") = BE - mHCO3/(1-Hct);
     Real BEe( unit = "meq/l")= BE + k*mHCO3/Hct;
     Real BE = 0;//time*40 - 20;
     Real mHCO3;

     Real SID;
     Real pCO2 = 40;
     Real Pi = 1.15;
     //  parameter Real alb = 4.4;
     Real alb = time*4.4 + 2;
     output Real pH = fullErythrocyte.pH;

    FiggeFencl3                         plasma_only(
       SID=normalPlasma.SID + BE,
       pCO2=pCO2,
       Pi=Pi,
       alb=alb) "Plasma model acc to FiggeFencl for comparison only" annotation (Placement(transformation(extent={{60,-20},{80,0}})));

     SAnomogram_formalization.SAoriginal normalBloodSA(
       Hct=Hct,
       BEox=BE,
       pCO2=pCO2) "SA original for comparison only"   annotation (Placement(transformation(extent={{60,20},{80,40}})));

    protected
     FiggeFenclNSID normalPlasma(
       pH0=7.4,
       pCO20=40,
       alb0=alb,
       Pi0=Pi)
       annotation (Placement(transformation(extent={{-80,40},{-60,60}})));

   //discrete Real dpic( start = 20);
   //discrete Real sTime( start = 0);
   equation
     plasma.pH = fullErythrocyte.pH;
     BEp = SID - normalPlasma.SID;
     annotation (experiment(
         StopTime=1,
         __Dymola_NumberOfIntervals=500,
         Tolerance=1e-003), Diagram(coordinateSystem(preserveAspectRatio=false,
             extent={{-100,-100},{100,100}}), graphics={
           Line(
             points={{-84,52},{-90,52},{-98,52},{-94,20},{-86,18}},
             color={28,108,200},
             arrow={Arrow.None,Arrow.Filled}),
           Line(
             points={{-54,16},{-30,2},{-26,-6}},
             color={28,108,200},
             arrow={Arrow.None,Arrow.Filled}),
           Line(
             points={{-54,-12},{-36,-14},{-26,-10}},
             color={28,108,200},
             arrow={Arrow.None,Arrow.Filled}),
           Text(
             extent={{-24,-14},{-4,0}},
             lineColor={28,108,200},
             fillColor={0,0,255},
             fillPattern=FillPattern.Solid,
             textString="pH"),
           Line(
             points={{10,84},{10,-74},{10,-80},{10,-82}},
             color={28,108,200},
             arrow={Arrow.None,Arrow.Filled})}));
   end test_combo_chng_Alb;

   model BErange
      import FullBloodAcidBase;
     Real logpCO2 = log10(pCO2);
     Real pCO2 = time*40 + 20;
     FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE_15(BE=-15,
          pCO2=pCO2)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
     FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE_10(BE=-10,
          pCO2=pCO2)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
     FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE_5(BE=-5,
          pCO2=pCO2)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
     FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE0(BE=-0,
          pCO2=pCO2)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
     FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE5(BE=5,
          pCO2=pCO2)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
     FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE10(BE=10,
          pCO2=pCO2)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
     FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE15(BE=15,
          pCO2=pCO2)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));

   /* 
// ALBUMIN CHARGE ACC TO FENCL OR OUR MODEL
createPlot(id=7, position={0, 0, 1350, 769}, y={"resultSetAtBE0.Normal.FF_plasma_only.AlbXMinus", "resultSetAtBE_15.Normal.FF_plasma_only.AlbXMinus",
 "resultSetAtBE15.Normal.FF_plasma_only.AlbXMinus", "resultSetAtBE_15.Normal.full_blood_combo.plasma.atch",
 "resultSetAtBE0.Normal.full_blood_combo.plasma.atch", "resultSetAtBE15.Normal.full_blood_combo.plasma.atch",
 "resultSetAtBE_15.Normal.FF_plasma_only.pH", "resultSetAtBE_15.Normal.full_blood_combo.pH",
 "resultSetAtBE0.Normal.full_blood_combo.pH", "resultSetAtBE0.Normal.FF_plasma_only.pH",
 "resultSetAtBE15.Normal.full_blood_combo.pH", "resultSetAtBE15.Normal.FF_plasma_only.pH",
 "resultSetAtBE_15.Normal.FF_plasma_only.atch", "resultSetAtBE0.Normal.FF_plasma_only.atch",
 "resultSetAtBE15.Normal.FF_plasma_only.atch"}, range={0.0, 1.0, -16.0, 9.0}, grid=true, filename="dsres.mat", colors={{238,46,47}, {28,108,200}, {0,140,72}, {180,56,148}, {0,0,0}, {162,29,33}, 
{244,125,35}, {102,44,145}, {28,108,200}, {238,46,47}, {0,140,72}, 
{180,56,148}, {0,0,0}, {162,29,33}, {244,125,35}}, patterns={LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Dash, LinePattern.Dash, LinePattern.Dash, LinePattern.Dash, 
LinePattern.Dash, LinePattern.Dash, LinePattern.Dash});
*/

   end BErange;

   model BErange_P
      import FullBloodAcidBase;
     Real logpCO2 = log10(pCO2);
     Real pCO2 = time*40 + 20;
     parameter Real Hb = 15;

     FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.ResultSetAtBE_P resultSetAtBE_10(
        BE=-10,
        pCO2=pCO2,
        Hb=Hb)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
     FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.ResultSetAtBE_P resultSetAtBE0(
        BE=-0,
        pCO2=pCO2,
        Hb=Hb)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
     FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.ResultSetAtBE_P resultSetAtBE10(
        BE=10,
        pCO2=pCO2,
        Hb=Hb)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));

       // PLOTS

       /* low, normal and high albumin for BE -10 0 10
createPlot(id=1, position={84, 99, 586, 418}, x="pCO2", y={"resultSetAtBE0.Normal.full_blood_combo.pH", "resultSetAtBE0.Normal.FF_plasma_only.pH",
 "resultSetAtBE0.LowP.full_blood_combo.pH", "resultSetAtBE0.HighP.FF_plasma_only.pH",
 "resultSetAtBE0.HighP.full_blood_combo.pH", "resultSetAtBE10.LowP.FF_plasma_only.pH",
 "resultSetAtBE10.Normal.full_blood_combo.plasma.pH", "resultSetAtBE10.Normal.full_blood_combo.fullErythrocyte.pH",
 "resultSetAtBE_10.Normal.full_blood_combo.plasma.pH", "resultSetAtBE_10.Normal.full_blood_combo.fullErythrocyte.pH",
 "resultSetAtBE_10.Normal.full_blood_combo.pH", "resultSetAtBE_10.HighP.full_blood_combo.pH",
 "resultSetAtBE_10.HighP.FF_plasma_only.pH", "resultSetAtBE10.Normal.FF_plasma_only.pH",
 "resultSetAtBE_10.Normal.FF_plasma_only.pH"}, range={20.0, 60.0, 7.0, 7.8}, grid=true, legend=false, filename="BErange_P.mat", logX=true, colors={{238,46,47}, {0,0,0}, {238,46,47}, {0,0,0}, {238,46,47}, {0,0,0}, {238,46,47}, 
{238,46,47}, {238,46,47}, {238,46,47}, {244,125,35}, {238,46,47}, {0,0,0}, 
{0,0,0}, {0,0,0}}, patterns={LinePattern.Solid, LinePattern.Solid, LinePattern.Dot, LinePattern.Dash, 
LinePattern.Dash, LinePattern.Dot, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Solid, LinePattern.Solid, LinePattern.Dash, LinePattern.Dash, 
LinePattern.Dash, LinePattern.Solid, LinePattern.Solid}, thicknesses={0.5, 0.25, 0.5, 0.0, 0.0, 0.25, 0.5, 0.5, 0.5, 0.5, 0.25, 0.0, 0.0, 0.25, 0.25});
 */

   end BErange_P;
 end Thrash;

  package SAnomogram_formalization

    partial model BaseAcidBase
     Real Hb(unit="g/dl") = Hct*33.34;
     input Real Hct(unit = "1") = 15/33.34;
     input Real BEox(displayUnit = "meq/l") = 0;
     input Real sO2(unit = "1") = 1;
     input Real pCO2(displayUnit = "mmHg") = 40;

    end BaseAcidBase;

    model Zander1995
      extends BaseAcidBase;
     //   Real Hb = Hct*33.34;
     // input Real Hct;
     // input Real BEox;
     // input Real sO2 = 1;
     // input Real pCO2;

      Real pH;
      Real BE = (1 - 0.0143*Hb)*(
                  (0.0304*pCO2*(10^(pH-6.1)) - 24.26)
                  + (9.5 + 1.63*Hb)*(pH - 7.4))
                  - 0.2*Hb*(1-sO2);
    equation
      BE = BEox + 0.2*(1-sO2)*Hb;
    end Zander1995;

    model Kofr2009
        extends BaseAcidBase;
     //   Real Hb = Hct*33.34;
     // input Real Hct;
     // input Real BEox;
     // input Real sO2 = 1;
     // input Real pCO2;

    Real a[:] = {996.35 - 10.35*Hb,
     35.16875 + 0.25875*Hb,
     -82.41 + 2.01*Hb,
     -5.27625 - (5.025e-2)*Hb,
     121 - Hb,
     2.625 + 0.025*Hb,
     -2.556 - 0.0944*Hb,
     13.87634 + 0.186653*Hb + 0.00534936*Hb^2,
     0.548 - 0.0274*Hb,
     0.274 - 0.0137*Hb};

    public
     Real BE = BEox + 0.2*(1-sO2)*Hb;
     Real Y = a[7] + sqrt(a[8] + a[9]*BE);
     Real pH = (a[1]*a[10] + a[2]*Y + (a[3]*a[10] + a[4]*Y)*log10(pCO2))/
                 (a[5]*a[10]  +a[6]*Y);

    end Kofr2009;

    model SAoriginal
          extends BaseAcidBase;
     //   Real Hb = Hct*33.34;
     // input Real Hct;
     // input Real BEox;
     // input Real sO2 = 1;
     // input Real pCO2;
    protected
      constant Real pco2BBCoef[:] = {2.1125e-009,  -640.9926e-009,     72.7649e-006,     -3.2862e-003,    -38.1749e-003,      8.2352e+000,    -97.0551e+000};
      constant Real pco2BECoef[:] = {8.3975e-009,   -513.9503e-009,      3.8105e-006,    231.6888e-006,    -46.5581e-003,    353.7105e-003,     39.9871e+000};
      constant Real pHBBCoef[:] = {40.8936e-012,    -13.0063e-009,      1.6780e-006,   -111.7919e-006,      4.0776e-003,    -67.8274e-003,      7.2888e+000};
      constant Real pHBECoef[:] = {131.3315e-009,      2.5027e-006,    175.6144e-006,     11.9273e-003,      7.4001e+000};
    public
      Real pCO2BB( start = 96) = pco2BBCoef[1]*BB^6 + pco2BBCoef[2]*BB^5 + pco2BBCoef[3]*BB^4 + pco2BBCoef[4]*BB^3 + pco2BBCoef[5]*BB^2  + pco2BBCoef[6]*BB + pco2BBCoef[7];
      Real pCO2BE( start = 40) = pco2BECoef[1]*BE^6 + pco2BECoef[2]*BE^5 + pco2BECoef[3]*BE^4 + pco2BECoef[4]*BE^3 + pco2BECoef[5]*BE^2  + pco2BECoef[6]*BE + pco2BECoef[7];
      Real pHBB( start = 7) = pHBBCoef[1]*BB^6 + pHBBCoef[2]*BB^5 + pHBBCoef[3]*BB^4 + pHBBCoef[4]*BB^3 + pHBBCoef[5]*BB^2  + pHBBCoef[6]*BB + pHBBCoef[7];
      Real pHBE( start = 7) = pHBECoef[1]*BE^4 + pHBECoef[2]*BE^3 + pHBECoef[3]*BE^2 + pHBECoef[4]*BE + pHBECoef[5];

    //   Real BB(start = 60) = BE+0.42*Hb+41.7;
    //   Real BE( start = -1)= BEox+0.2*(1-sO2)*Hb;
      Real BB = BE+0.42*Hb+41.7;
      Real BE = BEox+0.2*(1-sO2)*Hb;

      Real pH( start = 7.4) = (log10(pCO2) - log10(pCO2BB))*(pHBB - pHBE)/(log10(pCO2BB) - log10(pCO2BE))+ pHBB;
    equation

    end SAoriginal;

    model SAVanSlyke
      "The Van Slyke equation by Sigaard-Andersen, taken from http://www.siggaard-andersen.dk/OsaTextbook.htm"
      extends BaseAcidBase;
    //    parameter Real Hct = 15/33.34;
    //    Real BEox = 0;
    //    parameter Real sO2 = 1;
    //    parameter Real pCO2 = 40;

       Real ctH = -(1 - (1-rc)*fiEB)*((cHCO3 - cHCO30)
                   + buf*(pH - 7.4)) "concentration of titratable hydrogen ion";

       // Real rc = 0.57; // cHCO3e/cHCO3p from web pages
       Real rc = 0.57 - 0.28*dpH - 0.082*dpH^2; // said to be from siggaard-andersen, 1974: The acid-base status of the blood
       Real dpH = -(pH - 7.4);

       //parameter Real fiEB = 0.7;//15/33.34;// Hct;// ctHbB/cHBE
       Real fiEB = Hct;// ctHbB/cHBE
       Real ctHbE( unit = "mmol/l")= 21;
       Real cHCO30(unit = "mmol/l") = 24.5;
       Real buf = 2.3*ctHb + BetaP;
       Real ctHb = Hb*4 "SA is considering the monomer of hemoglobin tetramer"; //Hct*33.34; // recompute the hematocrit to mmol/l - BUT IS THAT COMPATIBLE? TODO!
       Real BetaP = 5.8 + 8.0 * (cAlb - 0.66);
       // Real BetaP = 7.7 + 8.0 * (cAlb - 0.66);
       Real cAlb( unit = "mmol/l")= Alb/6.6; //0.66;
       input Real Alb( unit = "g/dl") = 4.4;

       Real pH( start = 7.4)= 6.1 + log10(cHCO3/(0.230*pCO2*133/1000));
       Real BE;
       Real cHCO3(start = 24.5);// = 24.5;
    equation
      ctH = -BE;
      BEox = BE + 0.3*(1-sO2);

    end SAVanSlyke;

    model SAVanSlyke77
      "original Van Slyke equation according to SA 1977 DOI: 10.3109/00365517709098927"
      extends BaseAcidBase;
      Real a(unit = "mmol/l") "bicarbonate concentration in plasma/(mmol/l)";
      Real b(unit = "mmol/l") = Hb*0.6206
        "hemoglobin concentration in blood/(mmol/l)";
      Real c(unit = "1") = pH "pH of plasma at 37 degrees C";
      Real d(unit = "mEq/l") = BEox - 0.3*(1-sO2)
        "base excess concentration in blood/(mmol/l)";

      Real pH( start = 7.4)= 6.1 + log10(a/(0.230*pCO2*133/1000));
    equation
      a - 24.4 = - (2.3 * b + 7.7) * (c - 7.40) + d/(1 - 0.023 * b);
    end SAVanSlyke77;

    model SA_comparison_pCO2
      "pH dependent on varying pCO2 by different approximations. We take our implementation of SA nomogram as reference. Compare the object's pH"
      parameter Real Hct = 15/33.34;
      parameter Real BEox = 0;
      parameter Real sO2 = 1;
      output Real pHZander = zander1995.pH;
      output Real pHNomogram = sAoriginal.pH;
      output Real pHVanSlyke = sAVanSlyke.pH;
      input Real pCO2 = time*40 + 20;
      SAnomogram_formalization.Zander1995 zander1995(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
      SAnomogram_formalization.Kofr2009 kofr2009(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-20,40},{0,60}})));
      SAnomogram_formalization.SAoriginal sAoriginal(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{34,40},{54,60}})));
      SAnomogram_formalization.SAVanSlyke sAVanSlyke(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{80,40},{100,60}})));
      SAnomogram_formalization.SAVanSlyke77 sAVanSlyke77(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
      annotation (experiment(Tolerance=0.001), __Dymola_experimentSetupOutput,
        __Dymola_Commands(file="def.mos" "def"));
    end SA_comparison_pCO2;
  end SAnomogram_formalization;

  package Full_Blood

    model Full_blood_combo
      "Test Figge-fencl plasma and SA full hemoatocrite blood during variable pCO2"

      replaceable FiggeFencl3Detailed                         plasma(
        SID=SID,
        pCO2=pCO2,
        Pi=Pi,
        alb=alb)
        annotation (Placement(transformation(extent={{-20,20},{0,40}})));

      SAnomogram_formalization.SAoriginal fullErythrocyte(
        Hct=1,
        BEox=BEe,
        pCO2=pCO2)
        annotation (Placement(transformation(extent={{20,20},{40,40}})));
      constant Real fullHb = 33.34;
      input Real Hb = 15;
      input Real Hct = Hb/fullHb;

      Real BEp( unit = "meq/l") = BE - mHCO3/(1-Hct);
      Real BEe( unit = "meq/l")= BE + mHCO3/Hct;

      Real mHCO3;

      Real SID;
      input Real BE = 0;
      input Real pCO2 = 40;
      input Real Pi = 1.15;
      input Real alb = 4.4;
      output Real pH = fullErythrocyte.pH;

    protected
      FiggeFenclNSID normalPlasma(
        pH0=7.4,
        pCO20=40,
        alb0=alb,
        Pi0=Pi) "Computation of NSID"
        annotation (Placement(transformation(extent={{-40,60},{-20,80}})));
    equation
      plasma.pH = fullErythrocyte.pH;
      BEp = SID - normalPlasma.SID;
      annotation (experiment(
          StopTime=1,
          __Dymola_NumberOfIntervals=500,
          Tolerance=1e-003));
    end Full_blood_combo;

    package comparisson

      package Auxiliary

        model ResultSetAtBE
          parameter Real Pi = 1.15;
          parameter Real alb = 4.4;
          parameter Real BE = 0;
          parameter Real Hb = 15;
          output Real pHSABlood = normalBloodSA.pH;
          output Real pHComboBlood = Normal. full_blood_combo. pH;
          output Real pHFFplasma = Normal. FF_plasma_only. pH;
          output Real pHSAplasma = plasmaSA. pH;
          Real logpCO2 = log10(pCO2);
          input Real pCO2 = time*40 + 20;

          FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb LowAlb(  BE = BE, pCO2 = pCO2, alb = alb/2, Hb = Hb, Pi = Pi)
            "Low Alb"
          annotation (Placement(transformation(extent={{20,20},{40,40}})));

          FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb Normal(  BE = BE, pCO2 = pCO2, alb = alb, Hb = Hb, Pi = Pi)
            "Normal model set"
          annotation (Placement(transformation(extent={{20,20},{40,40}})));

          FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb HighAlb(  BE = BE, pCO2 = pCO2, alb = alb*2, Hb = Hb, Pi = Pi)
            "High Alb"
          annotation (Placement(transformation(extent={{20,20},{40,40}})));

        //   // BE 0
        //   full_blood_combo full_blood_combo1( BE = BE, pCO2 = pCO2, alb = 2, Hb = Hb) "Low Alb"
        //
        //   annotation (Placement(transformation(extent={{20,20},{40,40}})));
        //     full_blood_combo full_blood_combo2( BE = BE, pCO2 = pCO2, alb = alb, Hb = Hb) "Low Alb"
        //   annotation (Placement(transformation(extent={{20,20},{40,40}})));
        //
        //   full_blood_combo full_blood_combo3( BE = BE, pCO2 = pCO2, alb = 6.5, Hb = Hb) "Hogh Alb"
        //   annotation (Placement(transformation(extent={{20,20},{40,40}})));

        //   FiggeFencl3                         plasma_only(
        //     SID=normalPlasma.SID + BE,
        //     pCO2=pCO2,
        //     Pi=Pi,
        //     alb=alb) "Plasma model acc to FiggeFencl for comparison only" annotation (Placement(transformation(extent={{60,-20},{80,0}})));
        //
        //   FiggeFencl3                         plasma_low_alb(
        //     SID=lowAlbPlasma.SID + BE,
        //     pCO2=pCO2,
        //     Pi=Pi,
        //     alb=2.0) "Plasma model acc to FiggeFencl for comparison only with low albumin level" annotation (Placement(transformation(extent={{60,-20},{80,0}})));
        //
        //   FiggeFencl3                         plasma_high_alb(
        //     SID=highAlbPlasma.SID + BE,
        //     pCO2=pCO2,
        //     Pi=Pi,
        //     alb=6.5) "Plasma model acc to FiggeFencl for comparison only with high albumin level" annotation (Placement(transformation(extent={{60,-20},{80,0}})));

          SAnomogram_formalization.SAoriginal plasmaSA(
            Hb=0,
            BEox=BE,
            pCO2=pCO2) "SA original with Hct 0 for comparison with FF"   annotation (Placement(transformation(extent={{60,20},{80,40}})));

          SAnomogram_formalization.SAoriginal normalBloodSA(
            Hb = Hb,
            BEox=BE,
            pCO2=pCO2) "SA original for comparison only"   annotation (Placement(transformation(extent={{60,20},{80,40}})));

        // protected
        //   FiggeFenclNSID normalPlasma(
        //     pH0=7.4,
        //     pCO20=40,
        //     alb0=alb,
        //     Pi0=Pi)
        //     annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
        //   FiggeFenclNSID highAlbPlasma(
        //     pH0=7.4,
        //     pCO20=40,
        //     alb0=6.5,
        //     Pi0=Pi)
        //     annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
        //       FiggeFenclNSID lowAlbPlasma(
        //     pH0=7.4,
        //     pCO20=40,
        //     alb0=2.0,
        //     Pi0=Pi)
        //     annotation (Placement(transformation(extent={{-80,40},{-60,60}})));

        end ResultSetAtBE;

        model SetAtAlb

            Full_blood_combo full_blood_combo( BE = BE, pCO2 = pCO2, alb = alb, Hb = Hb, Pi=Pi)
          annotation (Placement(transformation(extent={{20,20},{40,40}})));

           FiggeFencl3Detailed                         FF_plasma_only(
            SID=normalPlasma.SID + BE,
            pCO2=pCO2,
            Pi=Pi,
            alb=alb) "Plasma model acc to FiggeFencl for comparison only" annotation (Placement(transformation(extent={{60,-20},{80,0}})));

            SAnomogram_formalization.SAoriginal normalBloodSA(
            Hb = Hb,
            BEox=BE,
            pCO2=pCO2) "SA original for comparison only"   annotation (Placement(transformation(extent={{60,20},{80,40}})));
            FullBloodAcidBase.SAnomogram_formalization.SAVanSlyke sAVanSlyke(  pCO2 = pCO2, BEox = BE, Hct = Hb/33.34, sO2 = 1, Alb= alb)
            annotation (Placement(transformation(extent={{-14,40},{6,60}})));
            Full_Blood.Wolf_full_blood wolf_full_blood(pCO2 = pCO2, BE = BE, AlbP = alb/AlbMW*1000*10/0.94);
            input Real BE, pCO2;
            input Real alb, Pi;
            input Real Hb;
        protected
            constant Real AlbMW = 66463;
            FiggeFenclNSID normalPlasma(
            pH0=7.4,
            pCO20=40,
            alb0=alb,
            Pi0=Pi)
            annotation (Placement(transformation(extent={{-80,40},{-60,60}})));

        end SetAtAlb;

        model ResultSetAtBE_P
          parameter Real Pi = 1.15;
          parameter Real alb = 4.4;
          parameter Real BE = 0;
          parameter Real Hb = 15;
          output Real pHSABlood = normalBloodSA.pH;
          output Real pHComboBlood = Normal. full_blood_combo. pH;
          output Real pHFFplasma = Normal. FF_plasma_only. pH;
          output Real pHSAplasma = plasmaSA. pH;
          Real logpCO2 = log10(pCO2);
          input Real pCO2 = time*40 + 20;

          FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb LowP(  BE = BE, pCO2 = pCO2, alb = alb, Hb = Hb, Pi = Pi/2)
            "Low Alb"
          annotation (Placement(transformation(extent={{20,20},{40,40}})));

          FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb Normal(  BE = BE, pCO2 = pCO2, alb = alb, Hb = Hb, Pi = Pi)
            "Normal model set"
          annotation (Placement(transformation(extent={{20,20},{40,40}})));

          FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb HighP(  BE = BE, pCO2 = pCO2, alb = alb, Hb = Hb, Pi = Pi*2)
            "High Alb"
          annotation (Placement(transformation(extent={{20,20},{40,40}})));

          SAnomogram_formalization.SAoriginal plasmaSA(
            Hb=0,
            BEox=BE,
            pCO2=pCO2) "SA original with Hct 0 for comparison with FF"   annotation (Placement(transformation(extent={{60,20},{80,40}})));

          SAnomogram_formalization.SAoriginal normalBloodSA(
            Hb = Hb,
            BEox=BE,
            pCO2=pCO2) "SA original for comparison only"   annotation (Placement(transformation(extent={{60,20},{80,40}})));

        end ResultSetAtBE_P;
      end Auxiliary;

    end comparisson;

    model Wolf_full_blood "Implementation of the Wolf acidbase model, reduced to Plasma and Erythrocyte"
      //default parameters
      //concentration in erythroctytes
      parameter Real NaE(unit="mmol/l") = 10/0.73;
      parameter Real KE(unit="mmol/l") = 99/0.73;
      parameter Real ClE(unit="mmol/l") = 53.8/0.73;
      parameter Real Hb(unit="mmol/l") = 5.3/0.73;
      parameter Real DPG(unit="mmol/l") = 4.4/0.73;
      parameter Real ATP(unit="mmol/l") = 1.8/0.73;
      parameter Real GSH(unit="mmol/l") = 2.2/0.73;
      parameter Real imE(unit="mmol/l") = 20.2/0.73;
      parameter Real PiE(unit="mmol/l") = 0.67/0.73;
      //concentrations in plasma nad intersticium
      parameter Real NaP(unit="mmol/l") = 140/0.94;
      parameter Real KP(unit="mmol/l") = 4.1/0.94;
      parameter Real CaP(unit="mmol/l") = 2.3/0.94;
      parameter Real MgP(unit="mmol/l") = 0.8/0.94;
      parameter Real ClP(unit="mmol/l") = 105/0.94;
      parameter Real PiP(unit="mmol/l") = 1.16/0.94;
      input Real AlbP(unit="mmol/l") = 0.65/0.94;
      parameter Real imP(unit="mmol/l") = 6/0.94;
      //charge on inpermeable solutes
      parameter Real ZimE = -9.2;
      parameter Real ZimP = -5.3;
      //volumes and CO2 preassure
      parameter Real Vblood(unit="l") = 5;
      input Real pCO2(unit="torr") = 20 + 40*time;
      //
      //derived parameters
      //water volumes in compartments
      parameter Real Vew0(unit="l") = 0.44 * 0.73 * Vblood;
      parameter Real Vpw0(unit="l") = (1 - 0.44) * 0.94 * Vblood;
      //masses in erythrocytes
      parameter Real m0NaE(unit="mmol") = NaE * Vew0;
      parameter Real m0KE(unit="mmol") = KE * Vew0;
      parameter Real m0ClE(unit="mmol") = ClE * Vew0;
      parameter Real m0Hb(unit="mmol") = Hb * Vew0;
      parameter Real m0DPG(unit="mmol") = DPG * Vew0;
      parameter Real m0ATP(unit="mmol") = ATP * Vew0;
      parameter Real m0GSH(unit="mmol") = GSH * Vew0;
      parameter Real m0imE(unit="mmol") = imE * Vew0;
      parameter Real m0PiE(unit="mmol") = PiE * Vew0;
      //masses in plasma
      parameter Real m0NaP(unit="mmol") = NaP * Vpw0;
      parameter Real m0KP(unit="mmol") = KP * Vpw0;
      parameter Real m0CaP(unit="mmol") = CaP * Vpw0;
      parameter Real m0MgP(unit="mmol") = MgP * Vpw0;
      parameter Real m0ClP(unit="mmol") = ClP * Vpw0;
      parameter Real m0PiP(unit="mmol") = PiP * Vpw0;
      Real m0AlbP(unit="mmol") = AlbP * Vpw0;
      parameter Real m0imP(unit="mmol") = imP * Vpw0;
      //overall masses of mobile ions
      Real MCl(unit="mmol");// = m0ClE + m0ClP;
      //
      //masses of mobile ions - 13 unknowns
      Real mClE(unit="mmol");
      Real mClP(unit="mmol",start = m0ClP);
      //
      //volumes of water - 3 unknowns
      Real Vew(unit="l",start = Vew0);
      Real Vpw(unit="l",start = Vpw0);
      //
      // concentrations of bicarbonates - 3 unknowns
      Real HCO3E(unit="mmol/l");
      Real HCO3P(unit="mmol/l");
      //
      //pH dependatn chargers
      Real ZPi = (-1) - 10 ^ (pHP - 6.87) / (1 + 10 ^ (pHP - 6.87));
      Real ZAlb = (-10.7) - 16 * (10 ^ (pHP - 7.42) / (1 + 10 ^ (pHP - 7.42)));
      Real ZHb = 15.6 - 23 * (10 ^ (pHE - 6.69) / (1 + 10 ^ (pHE - 6.69))) - 4 * (10 ^ (pHE - 7.89) / (1 + 10 ^ (pHE - 7.89))) + 1.5 * ((1 - 0.75) / 0.75);
      Real ZDPG = (-3) - 1 * (10 ^ (pHE - 7.56) / (1 + 10 ^ (pHE - 7.56))) - 1 * (10 ^ (pHE - 7.32) / (1 + 10 ^ (pHE - 7.32)));
      Real ZATP = (-3) - 1 * (10 ^ (pHE - 6.8) / (1 + 10 ^ (pHE - 6.8)));
      Real ZGSH = (-1) - 1 * (10 ^ (pHE - 8.54) / (1 + 10 ^ (pHE - 8.54))) - 1 * (10 ^ (pHE - 9.42) / (1 + 10 ^ (pHE - 9.42)));
      Real fiHb = 1 + 0.115 * C_Hb + 0.0256 * C_Hb ^ 2;
      //carbonates
      Real CO3E(unit="mmol/l") = HCO3E * 10 ^ (pHE - 10.2);
      Real CO3P(unit="mmol/l") = HCO3P * 10 ^ (pHP - 10.2);
      //concentrations
      Real C_NaE(unit="mmol/l",start = NaE) = m0NaE / Vew;
      Real C_KE(unit="mmol/l",start = KE) = m0KE / Vew;
      Real C_ClE(unit="mmol/l",start = ClE) = mClE / Vew;
      Real C_Hb(unit="mmol/l",start = Hb) = m0Hb / Vew;
      Real C_DPG(unit="mmol/l",start = DPG) = m0DPG / Vew;
      Real C_ATP(unit="mmol/l",start = ATP) = m0ATP / Vew;
      Real C_GSH(unit="mmol/l",start = GSH) = m0GSH / Vew;
      Real C_imE(unit="mmol/l",start = imE) = m0imE / Vew;
      Real C_PiE(unit="mmol/l",start = PiE) = m0PiE / Vew;
      //
      Real C_NaP(unit="mmol/l",start = NaP) = m0NaP / Vpw;
      Real C_KP(unit="mmol/l",start = KP) = m0KP / Vpw;
      Real C_CaP(unit="mmol/l",start = CaP) = m0CaP / Vpw;
      Real C_MgP(unit="mmol/l",start = MgP) = m0MgP / Vpw;
      Real C_ClP(unit="mmol/l",start = ClP) = mClP / Vpw;
      Real C_PiP(unit="mmol/l",start = PiP) = m0PiP / Vpw;
      Real C_AlbP(unit="mmol/l",start = AlbP) = m0AlbP / Vpw;
      Real C_imP(unit="mmol/l",start = imP) = m0imP / Vpw;
      //osmolality
      Real Oe = (0.93 * C_NaE + 0.93 * C_KE + 0.93 * C_ClE + 0.93 * C_PiE + fiHb * C_Hb + C_DPG + C_ATP + C_GSH + C_imE + HCO3E + CO3E);//*0.73;
      Real Op = (0.93 * C_NaP + 0.93 * C_KP + 0.93 * C_ClP + C_CaP + C_MgP + HCO3P + CO3P + 0.93 * C_PiP + C_AlbP + C_imP);// *0.94;
      //electric charge
      Real Qe = (C_NaE + C_KE - C_ClE - HCO3E - 2 * CO3E + ZHb * C_Hb + ZDPG * C_DPG + ZATP * C_ATP + ZGSH * C_GSH + ZimE)*Vew;
      Real Qp = (C_NaP + C_KP + 2 * C_CaP + 2 * C_MgP - C_ClP - HCO3P - 2 * CO3P + ZPi * C_PiP + ZAlb * C_AlbP + ZimP)*Vpw+X;
      parameter Real X(unit="mEq")=0;
      //
      input Real BE = 0;
     //Real SID;
      Real fH = (Vew + Vew / 0.73 * (1 - 0.73)) / (Vew + Vew / 0.73 * (1 - 0.73) + Vpw + Vpw / 0.94 * (1 - 0.94));
      Real rCl = C_ClE / C_ClP;
      Real rHCO3 = HCO3E / HCO3P;
      //addition of species
      parameter Real XCl(unit="mmol") = 0;

      Real HE(start=10^(-7.2));
      Real HP(start=10^(-7.37));
      //other unknowns
      //pH
      Real pHE(start = 7.22)=-log10(HE);
      Real pHP(start = 7.4)=-log10(HP);
    equation
      //mass conservation - 6 equations
      MCl = mClE + mClP - XCl;
      //volume conservation - 1 equation
      Vew0 + Vpw0 = Vew + Vpw;
      //
      //donnan equilibrium - 7 equations
      C_ClE / C_ClP = HP / HE;
      //
      //electroneutrality - 3 equations
      Qe = 0;
      Qp = 0;
        //osmotic equilibrium - 2 equations
      Op=Oe;
      //
      //carbonates and pH
      HCO3E = 0.026 * pCO2 * 10 ^ (pHE - 6.11);
      HCO3P = 0.0306 * pCO2 * 10 ^ (pHP - 6.11);
      //
      BE = (1 - 0.023 * 9) * (HCO3P - 24.4 + (2.3 * 9 + 7.7) * (pHP - 7.4));
      //SID = (1 - (1 - HCO3E / HCO3P) * fH * fB) * HCO3P + (1 - fH * fB) * (C_AlbP * (8 * pHP - 41) + C_PiP * (0.3 * pHP - 0.4)) + C_Hb * fB * (10.2 * pHP - 73.6) + C_DPG * fH * fB * (0.7 * pHP - 0.5);
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end Wolf_full_blood;
  end Full_Blood;

  package Figures

    model Figure1 "Comparison of SA nomogram formalisations"
      Real pCO2 = time*60 + 20;
      parameter Real Hb = 15;
      parameter Real Hct = Hb/33.34;
      SAnomogram_formalization.SA_comparison_pCO2 sA_comparison_pCO2_BE_10(BEox=
           -15, pCO2=pCO2, Hct = Hct)
        annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
      SAnomogram_formalization.SA_comparison_pCO2 sA_comparison_pCO2_BE0(BEox=0,
          pCO2=pCO2, Hct = Hct)
        annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
      SAnomogram_formalization.SA_comparison_pCO2 sA_comparison_pCO2_BE10(BEox=
            15, pCO2=pCO2, Hct = Hct)
        annotation (Placement(transformation(extent={{0,20},{20,40}})));
      FullBloodAcidBase.Figures.Figure2_4_BE_curve BE_curve;

      // Figure 1 Dymola Script
     /*
 createPlot(id=1, position={522, 188, 483, 300}, x="pCO2", y={"sA_comparison_pCO2_BE_10.sAoriginal.pH", "sA_comparison_pCO2_BE0.sAoriginal.pH",
 "sA_comparison_pCO2_BE10.sAoriginal.pH", "sA_comparison_pCO2_BE_10.zander1995.pH",
 "sA_comparison_pCO2_BE0.zander1995.pH", "sA_comparison_pCO2_BE10.zander1995.pH",
 "sA_comparison_pCO2_BE10.sAVanSlyke.pH", "sA_comparison_pCO2_BE0.sAVanSlyke.pH",
 "sA_comparison_pCO2_BE_10.sAVanSlyke.pH"}, range={20.0, 80.0, 7.0, 7.8}, erase=false, autoscale=false, grid=true, legend=false, filename="dsres.mat", logX=true, leftTitleType=0, bottomTitleType=0, colors={{0,140,72}, {0,140,72}, {0,140,72}, {0,0,0}, {0,0,0}, {0,0,0}, {28,108,200}, 
{28,108,200}, {28,108,200}}, patterns={LinePattern.Dash, LinePattern.Dash, LinePattern.Dash, LinePattern.Dot, 
LinePattern.Dot, LinePattern.Dot, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Solid}, thicknesses={0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25}, rightTitleType=0);
 
  */
    end Figure1;

    model Figure2_4_BE_curve "BE curve for graphs with SA nomogram"
      Real BE = time*60 - 30;
      output Real pCO2BE( start = 40) = pco2BECoef[1]*BE^6 + pco2BECoef[2]*BE^5 + pco2BECoef[3]*BE^4 + pco2BECoef[4]*BE^3 + pco2BECoef[5]*BE^2  + pco2BECoef[6]*BE + pco2BECoef[7];
      output Real pHBE( start = 7) = pHBECoef[1]*BE^4 + pHBECoef[2]*BE^3 + pHBECoef[3]*BE^2 + pHBECoef[4]*BE + pHBECoef[5];

    protected
      constant Real pco2BECoef[:] = {8.3975e-009,   -513.9503e-009,      3.8105e-006,    231.6888e-006,    -46.5581e-003,    353.7105e-003,     39.9871e+000};
      constant Real pHBECoef[:] = {131.3315e-009,      2.5027e-006,    175.6144e-006,     11.9273e-003,      7.4001e+000};

      // Dymola plot script

      /*
  createPlot(id=2, position={219, 4, 483, 300}, x="pCO2BE", y={"pHBE"}, range={20.0, 60.00000000000001, 7.1000000000000005, 7.800000000000001}, autoscale=false, grid=true, legend=false, filename="Figure1_4_BE_curve.mat", logX=true, leftTitleType=0, bottomTitleType=0, colors={{28,108,200}}, rightTitleType=0);
  */

    end Figure2_4_BE_curve;

    model figure3A "Verification of the Method during variable BE"

    //   Real BetaPlasma = der(plasma_only.pH)/der(BE);
    //   Real BetaEry = der(normalBloodSA.pH)/der(BE);
    //   Real BetaBlood = der(plasma.pH)/der(BE);
    /*
  FiggeFencl3Detailed                         plasma(
    SID=SID,
    pCO2=pCO2,
    Pi=Pi,
    alb=alb) "Plasma compartment (Hct = 0)"
    annotation (Placement(transformation(extent={{-80,8},{-60,28}})));
  */
      SAnomogram_formalization.SAoriginal plasma(
        Hct=0,
        BEox=BEp,
        pCO2=pCO2) "No erythrocyte blood with Hct = 0"
        annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

      SAnomogram_formalization.SAoriginal fullErythrocyte(
        Hct=1,
        BEox=BEe,
        pCO2=pCO2) "Full erythrocyte blood with Hct = 1"
        annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

      constant Real fullHb = 33.34;
      parameter Real Hb = 15;
      parameter Real Hct = Hb/fullHb;

      Real BEp( unit = "meq/l") = BE - mHCO3/(1-Hct);
      Real BEe( unit = "meq/l")= BE + mHCO3/Hct;
      //Real BE;
    //   Real BE = 0;
      Real BE = time*60 - 30;
      Real mHCO3;

    //  Real SID = 39;
     Real SID;
    //  Real pCO2 = time*40 + 20;
      Real pCO2 = 40;
      Real Pi = 1.15;
      parameter Real alb = 4.4;

     FiggeFencl3Detailed                         plasma_only(
        SID=normalPlasma.SID + BE,
        pCO2=pCO2,
        Pi=Pi,
        alb=alb) "Plasma model acc to FiggeFencl for comparison only" annotation (Placement(transformation(extent={{60,-20},{80,0}})));

      SAnomogram_formalization.SAoriginal normalBloodSA(
        Hct=Hct,
        BEox=BE,
        pCO2=pCO2) "SA original for comparison only"   annotation (Placement(transformation(extent={{60,20},{80,40}})));

    protected
      FiggeFenclNSID normalPlasma(
        pH0=7.4,
        pCO20=40,
        alb0=alb,
        Pi0=Pi)
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
    equation
      plasma.pH = fullErythrocyte.pH;
      BEp = SID - normalPlasma.SID;

    // Graph of varying BE
    /*
createPlot(id=21, position={0, 0, 483, 300}, x="BE", y={"plasma.pH", "normalBloodSA.pH"}, range={-30.0, 30.0, 6.6000000000000005, 7.800000000000001}, grid=true, legend=false, filename="figure2A.mat", leftTitleType=0, bottomTitleType=0, colors={{238,46,47}, {0,0,0}}, patterns={LinePattern.Solid, LinePattern.Dash}, thicknesses={0.5, 0.25}, rightTitleType=0);
*/
      annotation (experiment(StopTime=0.94),
                             Diagram(coordinateSystem(preserveAspectRatio=false,
              extent={{-100,-100},{100,100}}), graphics={
            Line(
              points={{-84,52},{-90,52},{-98,52},{-94,20},{-86,18}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Filled}),
            Line(
              points={{-54,16},{-30,2},{-26,-6}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Filled}),
            Line(
              points={{-54,-12},{-36,-14},{-26,-10}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Filled}),
            Text(
              extent={{-24,-14},{-4,0}},
              lineColor={28,108,200},
              fillColor={0,0,255},
              fillPattern=FillPattern.Solid,
              textString="pH"),
            Line(
              points={{10,84},{10,-74},{10,-80},{10,-82}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Filled})}),
        __Dymola_experimentSetupOutput);
    end figure3A;

    model figure3B "Verification of the Method during variable pCO2"

    //   Real BetaPlasma = der(plasma_only.pH)/der(BE);
    //   Real BetaEry = der(normalBloodSA.pH)/der(BE);
    //   Real BetaBlood = der(plasma.pH)/der(BE);
    /*
  FiggeFencl3Detailed                         plasma(
    SID=SID,
    pCO2=pCO2,
    Pi=Pi,
    alb=alb) "Plasma compartment (Hct = 0)"
    annotation (Placement(transformation(extent={{-80,8},{-60,28}})));
  */
      SAnomogram_formalization.SAoriginal plasma(
        Hct=0,
        BEox=BEp,
        pCO2=pCO2) "No erythrocyte blood with Hct = 0"
        annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

      SAnomogram_formalization.SAoriginal fullErythrocyte(
        Hct=1,
        BEox=BEe,
        pCO2=pCO2) "Full erythrocyte blood with Hct = 1"
        annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

      constant Real fullHb = 33.34;
      parameter Real Hb = 15;
      parameter Real Hct = Hb/fullHb;

      Real BEp( unit = "meq/l") = BE - mHCO3/(1-Hct);
      Real BEe( unit = "meq/l")= BE + mHCO3/Hct;
      //Real BE;
      Real BE = 0;
    //   Real BE = time*40 - 20;
      Real mHCO3;

    //  Real SID = 39;
     Real SID;
     Real pCO2 = time*40 + 20;
    //  Real pCO2 = 40;
      Real Pi = 1.15;
      parameter Real alb = 4.4;

     FiggeFencl3Detailed                         plasma_only(
        SID=normalPlasma.SID + BE,
        pCO2=pCO2,
        Pi=Pi,
        alb=alb) "Plasma model acc to FiggeFencl for comparison only" annotation (Placement(transformation(extent={{60,-20},{80,0}})));

      SAnomogram_formalization.SAoriginal normalBloodSA(
        Hct=Hct,
        BEox=BE,
        pCO2=pCO2) "SA original for comparison only"   annotation (Placement(transformation(extent={{60,20},{80,40}})));

    protected
      FiggeFenclNSID normalPlasma(
        pH0=7.4,
        pCO20=40,
        alb0=alb,
        Pi0=Pi)
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
    equation
      plasma.pH = fullErythrocyte.pH;
      BEp = SID - normalPlasma.SID;

      // grapf of varying pCO2
      /*
createPlot(id=22, position={95, 253, 483, 300}, x="pCO2", y={"plasma.pH", "normalBloodSA.pH"}, range={20.0, 60.0, 7.25, 7.65}, grid=true, legend=false, filename="figure2B.mat", leftTitleType=0, bottomTitleType=0, colors={{238,46,47}, {0,0,0}}, patterns={LinePattern.Solid, LinePattern.Dash}, thicknesses={1.0, 0.5}, rightTitleType=0);
*/
      annotation (experiment(
          StopTime=1,
          __Dymola_NumberOfIntervals=500,
          Tolerance=1e-003), Diagram(coordinateSystem(preserveAspectRatio=false,
              extent={{-100,-100},{100,100}}), graphics={
            Line(
              points={{-84,52},{-90,52},{-98,52},{-94,20},{-86,18}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Filled}),
            Line(
              points={{-54,16},{-30,2},{-26,-6}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Filled}),
            Line(
              points={{-54,-12},{-36,-14},{-26,-10}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Filled}),
            Text(
              extent={{-24,-14},{-4,0}},
              lineColor={28,108,200},
              fillColor={0,0,255},
              fillPattern=FillPattern.Solid,
              textString="pH"),
            Line(
              points={{10,84},{10,-74},{10,-80},{10,-82}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Filled})}));
    end figure3B;

    model Figure4
      Real logpCO2 = log10(pCO2);
      Real pCO2 = time*60 + 20;
      parameter Real Hb = 15;

      Full_Blood.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE_10(
        BE=-15,
        pCO2=pCO2,
        Hb=Hb)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
      Full_Blood.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE0(
        BE=-0,
        pCO2=pCO2,
        Hb=Hb)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
      Full_Blood.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE10(
        BE=15,
        pCO2=pCO2,
        Hb=Hb)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));

        // PLOTS
        /*
    fig 2 - Full blood combination compared to full blood SA, Figge of plasma and Wolf model
    createPlot(id=3, position={192, 255, 483, 300}, x="pCO2", y={"resultSetAtBE_10.Normal.FF_plasma_only.pH", "resultSetAtBE_10.Normal.full_blood_combo.pH",
 "resultSetAtBE_10.normalBloodSA.pH", "resultSetAtBE0.Normal.FF_plasma_only.pH",
 "resultSetAtBE0.Normal.full_blood_combo.pH", "resultSetAtBE0.normalBloodSA.pH",
 "resultSetAtBE10.Normal.FF_plasma_only.pH", "resultSetAtBE10.Normal.full_blood_combo.pH",
 "resultSetAtBE10.normalBloodSA.pH", "resultSetAtBE0.Normal.wolf_full_blood.pHP",
 "resultSetAtBE_10.Normal.wolf_full_blood.pHP", "resultSetAtBE10.Normal.wolf_full_blood.pHP"}, range={20.0, 80.0, 7.0, 7.8}, erase=false, autoscale=false, grid=true, legend=false, filename="dsres.mat", logX=true, leftTitleType=0, bottomTitleType=0, colors={{0,0,0}, {238,46,47}, {0,140,72}, {0,0,0}, {238,46,47}, {0,140,72}, {0,0,0}, 
{238,46,47}, {0,140,72}, {217,67,180}, {217,67,180}, {217,67,180}}, patterns={LinePattern.Solid, LinePattern.Solid, LinePattern.Dash, LinePattern.Solid, 
LinePattern.Solid, LinePattern.Dash, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Dash, LinePattern.Dot, LinePattern.Dot, LinePattern.Dot}, thicknesses={0.25, 0.5, 0.5, 0.25, 0.5, 0.5, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5}, rightTitleType=0);
    
    */

        /* FIGURE 3: Full blood combination compared to full blood SA, Figge of plasma and SA of Hct 0 (thus plasma as well)
    createPlot(id=3, position={70, 340, 483, 300}, x="pCO2", y={"resultSetAtBE_10.plasmaSA.pH", "resultSetAtBE_10.Normal.FF_plasma_only.pH", 
"resultSetAtBE_10.Normal.full_blood_combo.pH", "resultSetAtBE_10.normalBloodSA.pH",
 "resultSetAtBE0.plasmaSA.pH", "resultSetAtBE0.Normal.FF_plasma_only.pH", 
"resultSetAtBE0.Normal.full_blood_combo.pH", "resultSetAtBE0.normalBloodSA.pH",
 "resultSetAtBE10.plasmaSA.pH", "resultSetAtBE10.Normal.FF_plasma_only.pH", 
"resultSetAtBE10.Normal.full_blood_combo.pH", "resultSetAtBE10.normalBloodSA.pH"}, range={20.0, 80.0, 7.000000000000005, 7.799800000000000001}, autoscale=false, grid=true, legend=false, filename="BErange_alb_Hb5.mat", logX=true, leftTitleType=0, bottomTitleType=0, colors={{0,0,0}, {0,0,0}, {238,46,47}, {0,140,72}, {0,0,0}, {0,0,0}, {238,46,47}, 
{0,140,72}, {0,0,0}, {0,0,0}, {238,46,47}, {0,140,72}}, patterns={LinePattern.Dash, LinePattern.Solid, LinePattern.Solid, LinePattern.Dash, 
LinePattern.Dash, LinePattern.Solid, LinePattern.Solid, LinePattern.Dash, 
LinePattern.Dash, LinePattern.Solid, LinePattern.Solid, LinePattern.Dash}, thicknesses={0.25, 0.25, 0.5, 0.5, 0.25, 0.25, 0.5, 0.5, 0.25, 0.25, 0.5, 0.5}, rightTitleType=0);
*/

    /* FIGURE 4: low, normal and high albumin for BE -15 0 15
createPlot(id=4, position={0, 0, 483, 300}, x="pCO2", y={"resultSetAtBE0.Normal.full_blood_combo.pH", "resultSetAtBE0.Normal.sAVanSlyke.pH",
 "resultSetAtBE0.Normal.FF_plasma_only.pH", "resultSetAtBE0.LowAlb.sAVanSlyke.pH",
 "resultSetAtBE0.LowAlb.FF_plasma_only.pH", "resultSetAtBE0.LowAlb.full_blood_combo.pH",
 "resultSetAtBE0.HighAlb.sAVanSlyke.pH", "resultSetAtBE0.HighAlb.FF_plasma_only.pH",
 "resultSetAtBE0.HighAlb.full_blood_combo.pH", "resultSetAtBE10.LowAlb.full_blood_combo.pH",
 "resultSetAtBE10.LowAlb.FF_plasma_only.pH", "resultSetAtBE10.LowAlb.sAVanSlyke.pH",
 "resultSetAtBE10.Normal.full_blood_combo.plasma.pH", "resultSetAtBE10.Normal.full_blood_combo.fullErythrocyte.pH",
 "resultSetAtBE10.Normal.full_blood_combo.pH", "resultSetAtBE10.Normal.FF_plasma_only.pH",
 "resultSetAtBE10.Normal.sAVanSlyke.pH", "resultSetAtBE10.HighAlb.full_blood_combo.pH",
 "resultSetAtBE10.HighAlb.FF_plasma_only.pH", "resultSetAtBE10.HighAlb.sAVanSlyke.pH",
 "resultSetAtBE_10.LowAlb.FF_plasma_only.pH", "resultSetAtBE_10.LowAlb.sAVanSlyke.pH",
 "resultSetAtBE_10.Normal.full_blood_combo.plasma.pH", "resultSetAtBE_10.Normal.full_blood_combo.fullErythrocyte.pH",
 "resultSetAtBE_10.Normal.full_blood_combo.pH", "resultSetAtBE_10.Normal.FF_plasma_only.pH",
 "resultSetAtBE_10.Normal.sAVanSlyke.pH", "resultSetAtBE_10.HighAlb.full_blood_combo.pH",
 "resultSetAtBE_10.HighAlb.FF_plasma_only.pH", "resultSetAtBE_10.HighAlb.sAVanSlyke.pH",
 "resultSetAtBE_10.LowAlb.full_blood_combo.pH"}, range={20.0, 80.0, 7.00000005, 7.79998001}, autoscale=false, grid=true, legend=false, filename="BErange_alb_Hb5.mat", logX=true, leftTitleType=0, bottomTitleType=0, colors={{238,46,47}, {28,108,200}, {0,0,0}, {28,108,200}, {0,0,0}, {238,46,47}, 
{28,108,200}, {0,0,0}, {238,46,47}, {238,46,47}, {0,0,0}, {28,108,200}, 
{238,46,47}, {238,46,47}, {238,46,47}, {0,0,0}, {28,108,200}, {238,46,47}, 
{0,0,0}, {28,108,200}, {0,0,0}, {28,108,200}, {238,46,47}, {238,46,47}, 
{244,125,35}, {0,0,0}, {28,108,200}, {238,46,47}, {0,0,0}, {28,108,200}, 
{238,46,47}}, patterns={LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, LinePattern.Dot, 
LinePattern.Dot, LinePattern.Dot, LinePattern.Dash, LinePattern.Dash, 
LinePattern.Dash, LinePattern.Dot, LinePattern.Dot, LinePattern.Dot, 
LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Solid, LinePattern.Dash, LinePattern.Dash, LinePattern.Dash, 
LinePattern.Dot, LinePattern.Dot, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Dash, LinePattern.Solid, LinePattern.Solid, LinePattern.Dash, 
LinePattern.Dash, LinePattern.Dash, LinePattern.Dot}, thicknesses={0.5, 0.25, 0.25, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.5, 0.25, 0.25, 0.5, 0.5, 0.5,
 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 
 0.0, 0.5}, rightTitleType=0);
 
 /* FIGURE 4: low and high albumin range for BE -15, 0, 15
 createPlot(id=4, position={38, 451, 483, 300}, x="pCO2", y={"resultSetAtBE0.LowAlb.sAVanSlyke.pH", "resultSetAtBE0.LowAlb.FF_plasma_only.pH",
 "resultSetAtBE0.LowAlb.full_blood_combo.pH", "resultSetAtBE0.HighAlb.sAVanSlyke.pH",
 "resultSetAtBE0.HighAlb.FF_plasma_only.pH", "resultSetAtBE0.HighAlb.full_blood_combo.pH",
 "resultSetAtBE10.LowAlb.full_blood_combo.pH", "resultSetAtBE10.LowAlb.FF_plasma_only.pH",
 "resultSetAtBE10.LowAlb.sAVanSlyke.pH", "resultSetAtBE10.HighAlb.full_blood_combo.pH",
 "resultSetAtBE10.HighAlb.FF_plasma_only.pH", "resultSetAtBE10.HighAlb.sAVanSlyke.pH",
 "resultSetAtBE_10.LowAlb.FF_plasma_only.pH", "resultSetAtBE_10.HighAlb.full_blood_combo.pH",
 "resultSetAtBE_10.HighAlb.FF_plasma_only.pH", "resultSetAtBE_10.HighAlb.sAVanSlyke.pH",
 "resultSetAtBE_10.LowAlb.full_blood_combo.pH", "resultSetAtBE_10.LowAlb.sAVanSlyke.pH"}, range={20.0, 80.0, 7.000000005, 7.800000001}, autoscale=false, grid=true, legend=false, filename="dsres.mat", logX=true, leftTitleType=0, bottomTitleType=0, colors={{28,108,200}, {0,0,0}, {238,46,47}, {28,108,200}, {0,0,0}, {238,46,47}, 
{238,46,47}, {0,0,0}, {28,108,200}, {238,46,47}, {0,0,0}, {28,108,200}, 
{0,0,0}, {238,46,47}, {0,0,0}, {28,108,200}, {238,46,47}, {28,108,200}}, patterns={LinePattern.Dash, LinePattern.Dash, LinePattern.Dash, LinePattern.Solid, 
LinePattern.Solid, LinePattern.Solid, LinePattern.Dash, LinePattern.Dash, 
LinePattern.Dash, LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Dash, LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Dash, LinePattern.Dash}, thicknesses={0.25, 0.25, 0.5, 0.25, 0.25, 0.5, 0.5, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.5, 
0.25, 0.25, 0.5, 0.25}, rightTitleType=0);

*/
    end Figure4;

    model Figure5
      parameter Real Pi = 1.15;
      parameter Real alb = 4.4;
      parameter Real Hb = 15;
      Real BE = time*40 - 20;
      Real BEfixed = 0;
      Real pCO2 = time*60 + 20;
      Real pCO2fixed = 40;

      FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb fixedBE(  BE = BEfixed, pCO2 = pCO2, alb = alb, Hb = Hb, Pi = Pi)
        "Normal model set";
      FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb fixedpCO2(  BE = BE, pCO2 = pCO2fixed, alb = alb, Hb = Hb, Pi = Pi)
        "Normal model set";

        // A) VARYING BE, FIXED PCO2
        /*
    
    createPlot(id-51, position={96, 378, 483, 300}, x="BE", y={"fixedpCO2.full_blood_combo.plasma.barGraphHCO3Alb2", "fixedpCO2.full_blood_combo.plasma.barGraphHCO3Alb3",
 "fixedpCO2.full_blood_combo.plasma.barGraphHCO3Alb1", "fixedpCO2.FF_plasma_only.barGraphHCO3Alb1",
 "fixedpCO2.FF_plasma_only.barGraphHCO3Alb2", "fixedpCO2.FF_plasma_only.barGraphHCO3Alb3"}, range={-20.0, 20.0, 0.0, 60.0}, grid=true, legend=false, filename="Figure5.mat", leftTitleType=0, bottomTitleType=0, colors={{238,46,47}, {238,46,47}, {238,46,47}, {0,0,0}, {0,0,0}, {0,0,0}}, patterns={LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, LinePattern.Dash, 
LinePattern.Dash, LinePattern.Dash}, thicknesses={0.5, 0.5, 0.5, 0.5, 0.5, 0.5}, rightTitleType=0);
    
    */

        // B) VARYUING pCO2, FIXED BE
        /*
    
    createPlot(id=52, position={96, 378, 483, 300}, x="pCO2", y={"fixedBE.full_blood_combo.plasma.barGraphHCO3Alb1", "fixedBE.full_blood_combo.plasma.barGraphHCO3Alb2",
 "fixedBE.full_blood_combo.plasma.barGraphHCO3Alb3", "fixedBE.FF_plasma_only.barGraphHCO3Alb1",
 "fixedBE.FF_plasma_only.barGraphHCO3Alb2", "fixedBE.FF_plasma_only.barGraphHCO3Alb3"}, range={20.0, 60.0, 0.0, 45.0}, autoscale=false, grid=true, legend=false, filename="dsres.mat", leftTitleType=0, bottomTitleType=0, colors={{238,46,47}, {238,46,47}, {238,46,47}, {0,0,0}, {0,0,0}, {0,0,0}}, patterns={LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, LinePattern.Dash, 
LinePattern.Dash, LinePattern.Dash}, thicknesses={0.5, 0.5, 0.5, 0.5, 0.5, 0.5}, rightTitleType=0);    
    */
      annotation (Placement(transformation(extent={{20,20},{40,40}})));

    end Figure5;

    model Figure6

    //   full_blood_combo full_blood_combo1 annotation (Placement(transformation(extent={{-78,0},{-58,20}})));
    //   FiggeFencl3Detailed figgeFencl3Detailed annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
    //   SAnomogram_formalization.SAoriginal sAoriginal annotation (Placement(transformation(extent={{0,0},{20,20}})));
    //
      parameter Real pCO2 = 40;
    //  Real BE = 38.8*dilutionFactor - normalPlasma.SID;
      Real BE= ndp.SID*dilutionFactor - normalPlasma.SID;
      Real newSID = ndp.SID*dilutionFactor;
      // mnoram normal SID is 38.97
      parameter Real alb = 4.4;
      parameter Real Pi = 1.15;
      Real dilutionFactor = 1.25 - time;
      Real Hb = 15*dilutionFactor;
      Real HbPerCent = Hb/15*100;
      Full_Blood.comparisson.Auxiliary.SetAtAlb dilluted(
        BE=BE,
        pCO2=pCO2*dilutionFactor,
        alb=alb*dilutionFactor,
        Hb=Hb,
        Pi=Pi*dilutionFactor) "Dilluted"
        annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
      Full_Blood.comparisson.Auxiliary.SetAtAlb compensatedpCO2(
        BE=BE,
        pCO2=pCO2,
        alb=alb*dilutionFactor,
        Hb=Hb,
        Pi=Pi*dilutionFactor) "Dilluted"
        annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
      Full_Blood.comparisson.Auxiliary.SetAtAlb nondilluted(
        BE=BE,
        pCO2=pCO2,
        alb=alb,
        Hb=Hb,
        Pi=Pi) "Non-dilluted"
        annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));

        FiggeFenclNSID plasma(
          pH0=7.4,
          pCO20=pCO2*dilutionFactor,
          alb0=alb*dilutionFactor,
          Pi0=Pi*dilutionFactor);

    protected
        FiggeFenclNSID normalPlasma(
          pH0=7.4,
          pCO20=40,
          alb0=alb*dilutionFactor,
          Pi0=Pi*dilutionFactor);

        FiggeFenclNSID ndp(
          pH0=7.4,
          pCO20=40,
          alb0=alb,
          Pi0=Pi);

    public
      FiggeFencl3 figgeFencl3(
        SID=39.5*dilutionFactor,
        pCO2=pCO2*dilutionFactor,
        Pi=Pi*dilutionFactor,
        alb=alb*dilutionFactor) annotation (Placement(transformation(extent={{-80,0},{-60,20}})));

    // dilution graph
    /* 
createPlot(id=6, position={46, 289, 586, 421}, x="dilutionFactor", y={"compensatedpCO2.FF_plasma_only.pH", "dilluted.FF_plasma_only.pH", 
"dilluted.full_blood_combo.pH", "compensatedpCO2.full_blood_combo.pH"}, range={0.5, 1.25, 7.1000000000000005, 7.500000000000001}, grid=true, legend=false, filename="dsres.mat", bottomTitleType=2, bottomTitle="Dilution factor", colors={{0,0,0}, {0,0,0}, {238,46,47}, {238,46,47}}, patterns={LinePattern.Solid, LinePattern.Dash, LinePattern.Dash, LinePattern.Solid}, thicknesses={0.5, 0.5, 0.5, 0.5});
*/

    // dilution comparison to Zander

    end Figure6;

    model Figure7_9 "Change of pH while changing SID"
     // extends FiggeFencl3Base;

      FiggeFencl3Detailed figgeFencl3Base(
        SID=SID,
        pCO2=pCO2,
        Pi=Pi,
        alb=alb) annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
      Real SID = time*160 + 39 - 80;
      parameter Real pCO2 = 40;
      parameter Real Pi = 1.15;
      parameter Real alb = 4.4;
      output Real pH = figgeFencl3Base.pH;

      // FIGURE 7 - HOW IS THE TOTAL POSSIBLE CHARGE OF ALBUMIN built up

      /*
  
  createPlot(id=7, position={46, 289, 586, 421}, x="figgeFencl3Base.pH", y={"figgeFencl3Base.barGraphAlb2", "figgeFencl3Base.barGraphAlb4", 
"figgeFencl3Base.barGraphAlb1", "figgeFencl3Base.barGraphAlb3"}, range={4.0, 8.0, 0.0, 150.0}, autoscale=false, grid=true, legend=false, filename="dsres.mat", leftTitleType=0, bottomTitleType=0, colors={{0,0,0}, {0,0,0}, {28,108,200}, {238,46,47}}, thicknesses={1.0, 1.0, 0.5, 0.5}, rightTitleType=0);
  
  */

      // FIGURE 9 - HOW IS THE TOTAL CHARGE OF ALB AND PI DIFFERENT FROM LINEAR ESTIMATE
      /*
  
  // ALBUMIN  
  createPlot(id=91, position={55, 237, 300, 300}, x="figgeFencl3Base.pH", y={"figgeFencl3Base.AlbXMinus", "figgeFencl3Base.atch"}, range={6.0, 8.0, -18.0, -2.0}, autoscale=false, grid=true, legend=false, filename="dsres.mat", leftTitleType=0, bottomTitleType=0, colors={{0,0,0}, {238,46,47}}, patterns={LinePattern.Dash, LinePattern.Solid}, thicknesses={0.25, 0.5}, range2={0.05, 0.45}, rightTitleType=0);
  
  // PHOSPHATE
  createPlot(id=92, position={55, 237, 300, 300}, x="figgeFencl3Base.pH", y={"figgeFencl3Base.PXminus", "figgeFencl3Base.Pminus"}, range={6.0, 8.0, -2.3000000000000003, -1.6}, autoscale=false, grid=true, legend=false, filename="dsres.mat", leftTitleType=0, bottomTitleType=0, colors={{0,0,0}, {238,46,47}}, patterns={LinePattern.Dash, LinePattern.Solid}, thicknesses={0.25, 0.5}, rightTitleType=0);
    
  // TOGETHER
   createPlot(id=93, position={55, 237, 300, 300}, x="figgeFencl3Base.pH", y={"figgeFencl3Base.totalDiff"}, range={6.0, 8.0, -1.0, 0.40000000000000013}, autoscale=false, grid=true, legend=false, filename="dsres.mat", leftTitleType=0, bottomTitleType=0, colors={{0,0,0}}, thicknesses={0.5}, range2={0.84, 0.875}, rightTitleType=0);

  */

      annotation (experiment(Tolerance=1e-006), __Dymola_experimentSetupOutput);
    end Figure7_9;

    model Figure8

      Real BEf = 0;
      Real BEsidf = normalPlasma.SID - albPlasma.SID;

      parameter Real pCO2 = 40;
      parameter Real Hb = 15;
      parameter Real Pi = 1.15;
      parameter Real alb0 = 4.4;
      Real alb = alb0/2 + alb0*time "(0.5 - 1.5)*alb";

      FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb fixedBE(  BE = BEf, pCO2 = pCO2, alb = alb, Hb = Hb, Pi = Pi)
        "Normal model set"
      annotation (Placement(transformation(extent={{20,20},{40,40}})));

      FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb fixedSID(  BE = BEsidf, pCO2 = pCO2, alb = alb, Hb = Hb, Pi = Pi)
        "Normal model set"
      annotation (Placement(transformation(extent={{20,20},{40,40}})));

    protected
          FiggeFenclNSID albPlasma(
          pH0=7.4,
          pCO20=40,
          alb0=alb,
          Pi0=Pi);

          FiggeFenclNSID normalPlasma(
          pH0=7.4,
          pCO20=40,
          alb0=alb0,
          Pi0=Pi);

          // plot
          /*
createPlot(id=8, position={-39, 211, 865, 629}, x="alb", y={"fixedSID.FF_plasma_only.barGraphHCO3Alb1", "fixedSID.FF_plasma_only.barGraphHCO3Alb2",
 "fixedSID.FF_plasma_only.barGraphHCO3Alb3", "fixedBE.FF_plasma_only.barGraphHCO3Alb1",
 "fixedBE.FF_plasma_only.barGraphHCO3Alb2", "fixedBE.FF_plasma_only.barGraphHCO3Alb3",
 "fixedBE.FF_plasma_only.pH", "fixedSID.FF_plasma_only.pH"}, range={2.2, 6.6000000000000005, 0.0, 60.0}, grid=true, legend=false, filename="dsres.mat", leftTitleType=2, leftTitle="Charge [mEq/l]", bottomTitleType=2, bottomTitle="Albumin concentration [g/dl]", colors={{0,0,0}, {0,0,0}, {0,0,0}, {238,46,47}, {238,46,47}, {238,46,47}, {238,46,47}, 
{0,0,0}}, patterns={LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, LinePattern.Dot, 
LinePattern.Dot, LinePattern.Dot, LinePattern.Dash, LinePattern.Dash}, thicknesses={0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0}, range2={6.0, 7.6}, rightTitleType=2, rightTitle="pH", axes={1, 1, 1, 1, 1, 1, 2, 2});
*/
    end Figure8;
  end Figures;

  package Tests

    model testBE
      "Test combination of plasma and full hematocrit blood against original SA during varying BE"

      constant Real fullHb = 33.34;
      parameter Real Hct = 15/fullHb;
      Real BEp( unit = "meq/l") = BE - mHCO3/(1-Hct);
      Real BEe( unit = "meq/l")= BE + mHCO3/Hct;
      Real mHCO3;
      //Real pHFullBlood;

      Real BE;
      Real pCO2;
      SAnomogram_formalization.Zander1995 plasma(
        Hct=0,
        BEox=BEp,
        pCO2=pCO2)
        annotation (Placement(transformation(extent={{-40,18},{-20,38}})));
      SAnomogram_formalization.Zander1995 ery(
        Hct=1,
        BEox=BEe,
        pCO2=pCO2)
        annotation (Placement(transformation(extent={{-42,54},{-22,74}})));
      SAnomogram_formalization.SAoriginal SABlood(
        Hct=Hct,
        BEox=BE,
        pCO2=pCO2)
        annotation (Placement(transformation(extent={{40,40},{60,60}})));
    equation
      plasma.pH = ery.pH;
    //   plasma.pCO2 = ;
    //   ery.pCO2 = pCO2;

      // Hct = time + 1e-6;
      //Hct = 20/fullHb;
      BE = time*60 -30;
      pCO2 = 40;
    //   pCO2 = time * 40 + 20;
      annotation (experiment(
          StopTime=0.94,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-006), Diagram(coordinateSystem(preserveAspectRatio=false,
              extent={{-100,-100},{100,100}}), graphics={
            Line(points={{4,58},{6,62},{10,70},{20,70},{26,64},{26,56},{8,44},{10,30},
                  {24,28},{30,34},{30,40}}, color={28,108,200}),
            Ellipse(extent={{16,20},{24,12}}, lineColor={28,108,200}),
            Line(points={{-34,44},{-32,44},{-22,44},{-28,44},{-28,50},{-28,40}},
                color={28,108,200}),
            Line(points={{0,52},{2,52},{-8,52}}, color={28,108,200}),
            Line(points={{-8,38},{4,38}}, color={28,108,200})}));
    end testBE;

    model testpCO2
      "Test combination of plasma and full hematocrit blood against original SA during varying pCO2"

      SAnomogram_formalization.SAoriginal plasma(
        Hct=0,
        BEox=BEp,
        pCO2=pCO2)
        annotation (Placement(transformation(extent={{-20,20},{0,40}})));
      SAnomogram_formalization.SAoriginal ery(
        Hct=1,
        BEox=BEe,
        pCO2=pCO2)
        annotation (Placement(transformation(extent={{20,20},{40,40}})));
      SAnomogram_formalization.SAoriginal SABlood(
        Hct=Hct,
        BEox=BE,
        pCO2=pCO2)
        annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
      constant Real fullHb = 33.34;
      parameter Real Hct = 15/33.34;
      Real BEp( unit = "meq/l") = BE - mHCO3/(1-Hct);
      Real BEe( unit = "meq/l")= BE + mHCO3/Hct;
      Real mHCO3;
      //Real pHFullBlood;

      Real BE;
      Real pCO2;

    equation
      plasma.pH = ery.pH;
      BE = 0;//time*60 -30;
    //  pCO2 = time*20 + 20;
       pCO2 = time * 40 + 20;
      annotation (experiment(
          StopTime=1,
          __Dymola_NumberOfIntervals=500,
          Tolerance=1e-003));
    end testpCO2;

    model test_combo
      "Test Figge-fencl plasma and SA full hemoatocrite blood during variable pCO2 against "

    //   Real BetaPlasma = der(plasma_only.pH)/der(BE);
    //   Real BetaEry = der(normalBloodSA.pH)/der(BE);
    //   Real BetaBlood = der(plasma.pH)/der(BE);
    /*
  FiggeFencl3Detailed                         plasma(
    SID=SID,
    pCO2=pCO2,
    Pi=Pi,
    alb=alb) "Plasma compartment (Hct = 0)"
    annotation (Placement(transformation(extent={{-80,8},{-60,28}})));
  */
      parameter Real safe_Hct = 1;
      SAnomogram_formalization.SAoriginal plasma(
        Hct=safe_Hct,
        BEox=BEp,
        pCO2=pCO2) "No erythrocyte blood with Hct = 0"
        annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

      SAnomogram_formalization.SAoriginal fullErythrocyte(
        Hct=1,
        BEox=BEe,
        pCO2=pCO2) "Full erythrocyte blood with Hct = 1"
        annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));

      constant Real fullHb = 33.34;
      parameter Real Hb = 15;
      parameter Real Hct = Hb/fullHb;

      Real BEp( unit = "meq/l") = BE - mHCO3/(1-Hct);
      Real BEe( unit = "meq/l")= BE + mHCO3/Hct;
      //Real BE;
      Real BE = 0;
    //   Real BE = time*40 - 20;
      Real mHCO3;

    //  Real SID = 39;
     Real SID;
     Real pCO2 = time*40 + 20;
    //  Real pCO2 = 40;
      Real Pi = 1.15;
      parameter Real alb = 4.4;

     FiggeFencl3Detailed                         plasma_only(
        SID=normalPlasma.SID + BE,
        pCO2=pCO2,
        Pi=Pi,
        alb=alb) "Plasma model acc to FiggeFencl for comparison only" annotation (Placement(transformation(extent={{60,-20},{80,0}})));

      SAnomogram_formalization.SAoriginal normalBloodSA(
        Hct=Hct,
        BEox=BE,
        pCO2=pCO2) "SA original for comparison only"   annotation (Placement(transformation(extent={{60,20},{80,40}})));

    protected
      FiggeFenclNSID normalPlasma(
        pH0=7.4,
        pCO20=40,
        alb0=alb,
        Pi0=Pi)
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
    equation
      plasma.pH = fullErythrocyte.pH;
      BEp = SID - normalPlasma.SID;

      // grapf of varying pCO2
      /*
   createPlot(id=7, position={75, 70, 586, 421}, x="pCO2", y={"plasma.barGraphHCO3Alb1", "plasma.barGraphHCO3Alb2", "plasma.barGraphHCO3Alb3",
  "plasma_only.barGraphHCO3Alb1", "plasma_only.barGraphHCO3Alb2",
 "plasma_only.barGraphHCO3Alb3"}, range={20.0, 65.0, 0.0, 45.0}, grid=true, legend=false, filename="dsres.mat", colors={{238,46,47}, {238,46,47}, {238,46,47}, {0,0,0}, {0,0,0}, {0,0,0}}, patterns={LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, LinePattern.Dash,
 LinePattern.Dash, LinePattern.Dash}, thicknesses={0.5, 0.5, 0.5, 0.25, 0.25, 0.25});
*/
      annotation (experiment(
          StopTime=1,
          __Dymola_NumberOfIntervals=500,
          Tolerance=1e-003), Diagram(coordinateSystem(preserveAspectRatio=false,
              extent={{-100,-100},{100,100}}), graphics={
            Line(
              points={{-84,52},{-90,52},{-98,52},{-94,20},{-86,18}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Filled}),
            Line(
              points={{-54,16},{-30,2},{-26,-6}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Filled}),
            Line(
              points={{-54,-12},{-36,-14},{-26,-10}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Filled}),
            Text(
              extent={{-24,-14},{-4,0}},
              lineColor={28,108,200},
              fillColor={0,0,255},
              fillPattern=FillPattern.Solid,
              textString="pH"),
            Line(
              points={{10,84},{10,-74},{10,-80},{10,-82}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Filled})}));
    end test_combo;

    model SA_Figge_Fit
      Real pCO2 = time*40 + 20;
      parameter Real Kx( min = 0, max = 1e-3, nominal = 1e-7) = 0.00000000011 "K fixed";
      parameter Real Kxh = 0 "K dependent on H+";
      parameter Real Xi( min = -10, max = 10) = 1 "charge";

      Real error = (sA_comparison_pCO2_BE0.plasma_only.pH - sA_comparison_pCO2_BE0.sAoriginal.pH)^2;
    //  Real intError;
      SA_Figge_comparison_pCO2                    sA_comparison_pCO2_BE_10(BEox=
           -15, pCO2=pCO2, Hct = 0,
        Kx = Kx,
        Xi = Xi,
        Kxh = Kxh)
        annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
      SA_Figge_comparison_pCO2 sA_comparison_pCO2_BE0(BEox=0,
          pCO2=pCO2, Hct = 0,
        Kx = Kx,
        Xi = Xi,
        Kxh = Kxh)
        annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
      SA_Figge_comparison_pCO2 sA_comparison_pCO2_BE10(BEox=
            15, pCO2=pCO2, Hct = 0,
        Kx = Kx,
        Xi = Xi,
        Kxh = Kxh)
        annotation (Placement(transformation(extent={{0,20},{20,40}})));
      FullBloodAcidBase.Figures.Figure1_3_4_BE_curve BE_curve;

      // Figure 1 Dymola Script
     /*
createPlot(id=1, position={0, 0, 483, 300}, x="pCO2", y={"sA_comparison_pCO2_BE_10.sAoriginal.pH", "sA_comparison_pCO2_BE0.sAoriginal.pH",
 "sA_comparison_pCO2_BE10.sAoriginal.pH", "sA_comparison_pCO2_BE_10.zander1995.pH",
 "sA_comparison_pCO2_BE0.zander1995.pH", "sA_comparison_pCO2_BE10.zander1995.pH",
 "sA_comparison_pCO2_BE10.sAVanSlyke.pH", "sA_comparison_pCO2_BE0.sAVanSlyke.pH",
 "sA_comparison_pCO2_BE_10.sAVanSlyke.pH"}, range={20.0, 60.00000000000001, 7.1000000000000005, 7.800000000000001}, autoscale=false, grid=true, legend=false, filename="dsres.mat", logX=true, leftTitleType=0, bottomTitleType=0, colors={{238,46,47}, {238,46,47}, {238,46,47}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, 
{0,0,0}, {0,0,0}}, patterns={LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, LinePattern.Dot, 
LinePattern.Dot, LinePattern.Dot, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Solid}, thicknesses={1.0, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25}, rightTitleType=0);
*/

    // fig 5 comparisson to Figge
    /*
createPlot(id=5, position={95, 90, 586, 421}, x="pCO2", y={"sA_comparison_pCO2_BE0.zander1995.pH", "sA_comparison_pCO2_BE0.sAoriginal.pH",
 "sA_comparison_pCO2_BE0.plasma_only.pH", "sA_comparison_pCO2_BE_10.plasma_only.pH",
 "sA_comparison_pCO2_BE_10.zander1995.pH", "sA_comparison_pCO2_BE_10.sAoriginal.pH",
 "sA_comparison_pCO2_BE10.plasma_only.pH", "sA_comparison_pCO2_BE10.zander1995.pH",
 "sA_comparison_pCO2_BE10.sAoriginal.pH"}, range={20.0, 60.0, 7.1000000000000005, 7.800000000000001}, autoscale=false, grid=true, filename="SA_FIgge_Fit.mat", logX=true, colors={{0,140,72}, {238,46,47}, {0,0,0}, {0,0,0}, {0,140,72}, {238,46,47}, {0,0,0}, 
{0,140,72}, {238,46,47}}, patterns={LinePattern.Dot, LinePattern.Dash, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Dot, LinePattern.Dash, LinePattern.Solid, LinePattern.Dot, 
LinePattern.Dash}, thicknesses={0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5});
*/
    equation
    //  der(intError) = error;
    end SA_Figge_Fit;

    model SA_Figge_comparison_pCO2
      "pH dependent on varying pCO2 by different approximations. We take our implementation of SA nomogram as reference. Compare the object's pH"
      parameter Real Hct = 15/33.34;
      parameter Real BEox = 0;
      parameter Real sO2 = 1;
      output Real pHZander = zander1995.pH;
      output Real pHNomogram = sAoriginal.pH;
      output Real pHVanSlyke = sAVanSlyke.pH;
      input Real pCO2 = time*40 + 20;
      parameter Real Pi = 1.15;
      parameter Real alb = 4.4;
      parameter Real Kx;
      parameter Real Xi;
      parameter Real Kxh;
      FiggeFencl3Extended                         plasma_only(
        SID=normalPlasma.SID + BEox,
        pCO2=pCO2,
        Pi=Pi,
        alb=alb,
        Kx = Kx,
        Xi = Xi,
        Kxh = Kxh) "Plasma model acc to FiggeFencl for comparison only"
                                                               annotation (Placement(transformation(extent={{60,-20},{80,0}})));

      SAnomogram_formalization.Zander1995 zander1995(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
      SAnomogram_formalization.Kofr2009 kofr2009(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-20,40},{0,60}})));
      SAnomogram_formalization.SAoriginal sAoriginal(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{34,40},{54,60}})));
      SAnomogram_formalization.SAVanSlyke sAVanSlyke(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{80,40},{100,60}})));
      SAnomogram_formalization.SAVanSlyke77 sAVanSlyke77(
        pCO2=pCO2,
        BEox=BEox,
        Hct=Hct,
        sO2=sO2)
        annotation (Placement(transformation(extent={{-60,40},{-40,60}})));

    protected
      FiggeFenclNSID normalPlasma(
        pH0=7.4,
        pCO20=40,
        alb0=alb,
        Pi0=Pi)
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
      annotation (experiment(Tolerance=0.001), __Dymola_experimentSetupOutput,
        __Dymola_Commands(file="def.mos" "def"));
    end SA_Figge_comparison_pCO2;

    model Test_Combo_Wolf
      Real pCO2 = 20 + time*60;
      Real BE = 0;
      parameter Real a = 0.65/0.94;
      Full_Blood.Full_blood_combo full_blood_combo( pCO2 = pCO2, BE = BE)
        annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
      Full_Blood.Wolf_full_blood wolf_full_blood(pCO2 = pCO2, BE = BE, AlbP = a)
        annotation (Placement(transformation(extent={{0,0},{20,20}})));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Test_Combo_Wolf;

    model Test_Combo_Wolf_at_BE
      Real pCO2 = 20 + time*60;
      Test_Combo_Wolf test_Combo_Wolf_15( pCO2 = pCO2, BE = -15)
        annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
      Test_Combo_Wolf test_Combo_Wolf0( pCO2 = pCO2, BE = 0)
        annotation (Placement(transformation(extent={{-20,0},{0,20}})));
      Test_Combo_Wolf test_Combo_Wolf15( pCO2 = pCO2, BE = 15)
        annotation (Placement(transformation(extent={{20,0},{40,20}})));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));

    /*        
// Dymola plot command

createPlot(id=1, position={15, 10, 584, 420}, x="pCO2", y={"test_Combo_Wolf_15.full_blood_combo.pH", "test_Combo_Wolf_15.wolf_full_blood.pHP", "test_Combo_Wolf0.full_blood_combo.pH", "test_Combo_Wolf0.wolf_full_blood.pHP", "test_Combo_Wolf15.wolf_full_blood.pHP", "test_Combo_Wolf15.full_blood_combo.pH"}, range={20.0, 80.0, 7.0, 7.8}, autoscale=false, grid=true, legend=false, filename="dsres.mat", logX=true, leftTitleType=2, leftTitle="pH", colors={{238,46,47}, {0,0,0}, {238,46,47}, {0,0,0}, {0,0,0}, {238,46,47}}, patterns={LinePattern.Solid, LinePattern.Dash, LinePattern.Solid, LinePattern.Dash, LinePattern.Dash, LinePattern.Solid}, thicknesses={0.5, 0.5, 0.5, 0.5, 0.5, 0.5});
*/
    end Test_Combo_Wolf_at_BE;
  end Tests;

  package AlbuminBorderFlux
    model AlbuminBalance
      import FullBloodAcidBase;

    Physiolibrary.Chemical.Components.Diffusion UT_Capillary(Conductance(
          displayUnit="l/day") = 3.9351851851852e-09)
      annotation (Placement(transformation(extent={{-2,86},{10,98}})));
    Physiolibrary.Chemical.Components.Diffusion MT_Capillary(Conductance(
          displayUnit="l/day") = 7.4074074074074e-09)
      annotation (Placement(transformation(extent={{-2,54},{10,66}})));
    Physiolibrary.Chemical.Components.Diffusion LT_Capillary(Conductance(
          displayUnit="l/day") = 1.1805555555556e-08)
      annotation (Placement(transformation(extent={{0,28},{12,40}})));
    Physiolibrary.Chemical.Sources.UnlimitedSolutePump Transfusion(
          useSoluteFlowInput=false, SoluteFlow=0)
      annotation (Placement(transformation(extent={{20,-22},{0,-2}})));
    Physiolibrary.Chemical.Components.Stream UT_Lymph(useSolutionFlowInput=
            false, SolutionFlow=5.5333333333333e-09)
      annotation (Placement(transformation(extent={{10,82},{0,72}})));
    Physiolibrary.Chemical.Components.Stream MT_Lymph(useSolutionFlowInput=
            false, SolutionFlow=1.315e-08)
      annotation (Placement(transformation(extent={{10,50},{0,40}})));
    Physiolibrary.Chemical.Components.Stream LT_Lymph(useSolutionFlowInput=
            false, SolutionFlow=1.5933333333333e-08)
      annotation (Placement(transformation(extent={{10,24},{0,14}})));
      Physiolibrary.Chemical.Components.Substance plasma(
        stateName="PlasmaProtein.Mass",
        useNormalizedVolume=false,
        solute_start=0.00406)
        annotation (Placement(transformation(extent={{-76,40},{-56,60}})));
      Physiolibrary.Chemical.Components.Substance UpperTorso(
        stateName="UT_InterstitialProtein.Mass",
        useNormalizedVolume=false,
        solute_start=0.00122)
        annotation (Placement(transformation(extent={{78,72},{58,92}})));
      Physiolibrary.Chemical.Components.Substance MiddleTorso(
        stateName="MT_InterstitialProtein.Mass",
        useNormalizedVolume=false,
      solute_start=0.00299)
        annotation (Placement(transformation(extent={{78,42},{58,62}})));
      Physiolibrary.Chemical.Components.Substance LowerTorso(
        stateName="LT_InterstitialProtein.Mass",
        useNormalizedVolume=false,
      solute_start=0.0018)
        annotation (Placement(transformation(extent={{78,14},{58,34}})));
      AlbuminSynthesis                               synthesis(SynthesisBasic=
            1.6666666666667e-07, UseSythesisFactorInput=false)
        annotation (Placement(transformation(extent={{-80,-66},{-60,-46}})));
      Degradation                                      degradation(
          DegradationBasic=1.6666666666667e-07, UseDegradationFactorInput=false)
        annotation (Placement(transformation(extent={{20,-54},{40,-34}})));
    Physiolibrary.Chemical.Components.Diffusion GlomerulusProtein_Perm(
        Conductance=(0)*(1e-6)/60)
      annotation (Placement(transformation(extent={{0,-8},{20,12}})));
      Physiolibrary.Chemical.Components.Substance Bladder(
        stateName="BladderProtein.Mass",
        useNormalizedVolume=false,
      solute_start=1e-15)
        annotation (Placement(transformation(extent={{78,-8},{58,12}})));
    Physiolibrary.Chemical.Sensors.ConcentrationMeasure concentrationMeasure1
      annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=180,
          origin={44,-24})));
      ProteinDivision proteinDivision
        annotation (Placement(transformation(extent={{52,-38},{60,-30}})));
      Physiolibrary.Types.Constants.VolumeConst volume(k=0.006063)
        annotation (Placement(transformation(extent={{98,52},{90,60}})));
      Physiolibrary.Types.Constants.VolumeConst volume1(k=0.00185)
        annotation (Placement(transformation(extent={{98,82},{90,90}})));
      Physiolibrary.Types.Constants.VolumeConst volume2(k=0.00247)
        annotation (Placement(transformation(extent={{98,24},{90,32}})));
      Physiolibrary.Types.Constants.VolumeConst volume3(k=0.0003)
        annotation (Placement(transformation(extent={{98,2},{90,10}})));
      Physiolibrary.Types.Constants.VolumeConst volume4(k=0.002807) annotation (
         Placement(transformation(
            extent={{-4,-4},{4,4}},
            rotation=0,
            origin={-74,68})));
      Physiolibrary.Chemical.Sensors.MolarFlowMeasure molarFlowMeasure
        annotation (Placement(transformation(extent={{-54,-66},{-34,-46}})));
      Physiolibrary.Chemical.Sensors.MolarFlowMeasure molarFlowMeasure1
        annotation (Placement(transformation(extent={{-22,-66},{-2,-46}})));
      Modelica.Blocks.Math.Add add(k1=-1, k2=+1)
        annotation (Placement(transformation(extent={{0,-86},{10,-76}})));
      ProteinDivision proteinDivision1
        annotation (Placement(transformation(extent={{20,-86},{30,-76}})));
      ProteinCharge proteinCharge(UseSimpleAlbuminCharge=true)
        annotation (Placement(transformation(extent={{40,-96},{60,-76}})));
      AcidBaseBuffers acidBaseBuffers(UseConstantAlb=true, plasmaVol(displayUnit="l"))
        annotation (Placement(transformation(extent={{80,-96},{100,-76}})));
      Physiolibrary.Chemical.Components.Clearance degradation1(
          useSolutionFlowInput=true)
        annotation (Placement(transformation(extent={{20,-70},{40,-50}})));
      FullBloodAcidBase.AlbuminBorderFlux.pulse pulse(
        normal=0,
        startTime(displayUnit="h") = 172800,
        length(displayUnit="h") = 86400,
        dose=16e-7)
        annotation (Placement(transformation(extent={{68,-64},{48,-44}})));
    equation
      connect(UT_Capillary.q_out,UpperTorso. q_out) annotation (Line(
          points={{10,92},{18,92},{18,82},{68,82}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,UT_Capillary. q_in) annotation (Line(
          points={{-66,50},{-26,50},{-26,92},{-2,92}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,UT_Lymph. q_out) annotation (Line(
          points={{-66,50},{-26,50},{-26,76},{0,76},{0,77}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,MT_Capillary. q_in) annotation (Line(
          points={{-66,50},{-26,50},{-26,58},{-14,58},{-14,60},{-2,60}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(MT_Capillary.q_out,MiddleTorso. q_out) annotation (Line(
          points={{10,60},{16,60},{16,52},{68,52}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(MT_Lymph.q_in,MiddleTorso. q_out) annotation (Line(
          points={{10,45},{16,45},{16,52},{68,52}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,MT_Lymph. q_out) annotation (Line(
          points={{-66,50},{-26,50},{-26,44},{0,44},{0,45}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,LT_Capillary. q_in) annotation (Line(
          points={{-66,50},{-26,50},{-26,34},{0,34}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,LT_Lymph. q_out) annotation (Line(
          points={{-66,50},{-26,50},{-26,18},{0,18},{0,19}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,GlomerulusProtein_Perm. q_in) annotation (Line(
          points={{-66,50},{-26,50},{-26,2},{0,2}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,Transfusion. q_out) annotation (Line(
          points={{-66,50},{-26,50},{-26,-12},{0,-12}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(LT_Capillary.q_out,LowerTorso. q_out) annotation (Line(
          points={{12,34},{18,34},{18,24},{68,24}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(LT_Lymph.q_in,LowerTorso. q_out) annotation (Line(
          points={{10,19},{18,19},{18,24},{68,24}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(GlomerulusProtein_Perm.q_out,Bladder. q_out) annotation (Line(
          points={{20,2},{68,2}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,concentrationMeasure1. q_in) annotation (Line(
          points={{-66,50},{-26,50},{-26,-24},{44,-24}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
    connect(concentrationMeasure1.concentration,proteinDivision. totalProteins)
      annotation (Line(
        points={{44,-32},{44,-34},{52,-34}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(UT_Lymph.q_in,UpperTorso. q_out) annotation (Line(
        points={{10,77},{18,77},{18,82},{68,82}},
        color={107,45,134},
        thickness=1,
        smooth=Smooth.None));
      connect(MiddleTorso.solutionVolume, volume.y)
        annotation (Line(points={{72,56},{72,56},{89,56}}, color={0,0,127}));
      connect(UpperTorso.solutionVolume, volume1.y)
        annotation (Line(points={{72,86},{89,86}}, color={0,0,127}));
      connect(LowerTorso.solutionVolume, volume2.y)
        annotation (Line(points={{72,28},{89,28}}, color={0,0,127}));
      connect(Bladder.solutionVolume, volume3.y) annotation (Line(points={{72,6},{
              72,6},{89,6}},           color={0,0,127}));
      connect(plasma.solutionVolume, volume4.y) annotation (Line(points={{-70,54},
              {-70,54},{-70,68},{-69,68}},     color={0,0,127}));
      connect(synthesis.q_out, molarFlowMeasure.q_in) annotation (Line(
          points={{-60,-56},{-54,-56}},
          color={107,45,134},
          thickness=1));
      connect(molarFlowMeasure.q_out, plasma.q_out) annotation (Line(
          points={{-34,-56},{-26,-56},{-26,50},{-66,50}},
          color={107,45,134},
          thickness=1));
      connect(add.u1, molarFlowMeasure1.molarFlowRate) annotation (Line(points={{-1,-78},
              {-12,-78},{-12,-64}},           color={0,0,127}));
      connect(add.u2, molarFlowMeasure.molarFlowRate) annotation (Line(points={{-1,-84},
              {-2,-84},{-14,-84},{-44,-84},{-44,-64}},color={0,0,127}));
      connect(add.y, proteinDivision1.totalProteins) annotation (Line(points={{10.5,
              -81},{10.5,-81},{20,-81}},      color={0,0,127}));
      connect(proteinCharge.port_a, acidBaseBuffers.port_a) annotation (Line(
          points={{59,-86},{59,-86},{81,-86}},
          color={107,45,134},
          thickness=1));
      connect(acidBaseBuffers.pH, proteinCharge.pH) annotation (Line(points={{
              81,-95},{81,-100},{36,-100},{36,-95},{40,-95}}, color={0,0,127}));
      connect(molarFlowMeasure1.q_in, plasma.q_out) annotation (Line(
          points={{-22,-56},{-26,-56},{-26,50},{-66,50}},
          color={107,45,134},
          thickness=1));
      connect(proteinDivision1.albumin, proteinCharge.AlbuminDifferenceMolarFlow)
        annotation (Line(points={{30,-78},{40,-78},{40,-77}}, color={0,0,127}));
      connect(proteinDivision.albumin, acidBaseBuffers.albuminConcentration)
        annotation (Line(points={{60,-31.6},{70,-31.6},{70,-32},{80,-32},{80,
              -77}}, color={0,0,127}));
      connect(molarFlowMeasure1.q_out, degradation.q_in) annotation (Line(
          points={{-2,-56},{10,-56},{10,-44},{20,-44}},
          color={107,45,134},
          thickness=1));
      connect(molarFlowMeasure1.q_out, degradation1.q_in) annotation (Line(
          points={{-2,-56},{10,-56},{10,-60},{20,-60}},
          color={107,45,134},
          thickness=1));
      connect(degradation1.solutionFlow, pulse.y) annotation (Line(points={{30,
              -53},{40,-53},{40,-54},{49,-54}}, color={0,0,127}));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}})), experiment(StopTime=1.728e+006,
            __Dymola_NumberOfIntervals=5000));
    end AlbuminBalance;

    model AlbuminSynthesis
    //  parameter Physiolibrary.Types.MassFlowRate  SynthesisBasic "10 mg/min";
      parameter Physiolibrary.Types.MolarFlowRate SynthesisBasic = 2.75753e-09
        "10 mg/min";
      parameter Real[:,3] data =  {{ 20.0,  3.0,  0.0}, { 28.0,  1.0,  -0.2}, { 40.0,  0.0,  0.0}}
        "COPEffect";
    Physiolibrary.Blocks.Interpolation.Curve c(
      x=data[:, 1],
      y=data[:, 2],
      slope=data[:, 3],
      Xscale=101325/760);

    Physiolibrary.Chemical.Interfaces.ChemicalPort_b q_out annotation (extent=[
          -10,-110; 10,-90], Placement(transformation(extent={{90,-10},{110,10}})));

      Physiolibrary.Types.Pressure COP;
    //  Physiolibrary.Types.AmountOfSubstance  synthetizedAmount(start=0);
    //  Physiolibrary.Types.Mass  synthetizedMass(start=0);
    //protected
    //  constant Physiolibrary.Types.Time sec=1;
    //  constant Physiolibrary.Types.Volume ghostPlasmaVol=3.02e-3
    //    "Strange dependence derived from original HumMod";
      Modelica.Blocks.Interfaces.RealInput SynthesisFactor = SyntFact if UseSythesisFactorInput
        annotation (Placement(transformation(extent={{-100,60},{-60,100}}),
            iconTransformation(extent={{-100,60},{-60,100}})));

    parameter Real SynthesisFactorParam = 1
     annotation (Dialog(enable=not UseSythesisFactorInput));
    parameter Boolean UseSythesisFactorInput = false
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="External inputs/outputs"));
    protected
      Real SyntFact;
    equation
      if not UseSythesisFactorInput then
        SyntFact = SynthesisFactorParam;
      end if;
      COP =  q_out.conc * Modelica.Constants.R * 310.15;
      c.u=COP;
      q_out.q = -SynthesisBasic * c.val*SyntFact;

    //TODO: state
    //der(synthetizedAmount) = -q_out.q;
    //  ProteinsMass2AmountOfSubstance(synthetizedMass,ghostPlasmaVol) = synthetizedAmount;
     annotation (
        defaultComponentName="synthesis",
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                            graphics={Rectangle(
              extent={{-100,-50},{100,50}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid), Text(
              extent={{-100,-50},{90,50}},
              lineColor={0,0,255},
              textString="%name")}),  Diagram(coordinateSystem(preserveAspectRatio=true,
                       extent={{-100,-100},{100,100}}), graphics),
        Documentation(revisions="<html>
<p><i>2009-2010</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics));
    end AlbuminSynthesis;

    model Degradation
    //  parameter Physiolibrary.Types.MassFlowRate  DegradationBasic "10 mg/min";
    //  parameter Real[:,3] data =  {{ 0.00,  0.0,  0.0}, { 0.07,  1.0,  40.0}, { 0.09,  6.0,  0.0}}
    //    "ProteinEffect";
       parameter Physiolibrary.Types.MolarFlowRate DegradationBasic = 2.75753e-09
        "10 mg/min";
       parameter Real[:,3] data =  {{ 0.00,  0.0,  0.0}, { 1.45,  1.0,  1.59}, { 1.97,  6.0,  0.0}}
        "ProteinEffect";

    Physiolibrary.Blocks.Interpolation.Curve c(
      x=data[:, 1],
      y=data[:, 2],
      slope=data[:, 3]);
    Physiolibrary.Chemical.Interfaces.ChemicalPort_a q_in annotation (Placement(
          transformation(extent={{-100,0},{-60,40}}), iconTransformation(extent=
             {{-110,-10},{-90,10}})));

    //  Physiolibrary.Types.AmountOfSubstance  degradedAmount(start=0);
    //  Physiolibrary.Types.Mass  degradedMass(start=0);
    //protected
    //  constant Physiolibrary.Types.Time sec=1;
    //  constant Physiolibrary.Types.Volume ghostPlasmaVol=3.02e-3
    //    "Strange dependence derived from original HumMod";
      Modelica.Blocks.Interfaces.RealInput DegradationFactor = DegFact if UseDegradationFactorInput
        annotation (Placement(transformation(extent={{-100,60},{-60,100}}),
            iconTransformation(extent={{-100,60},{-60,100}})));

    parameter Real DegradationFactorParam = 1
     annotation (Dialog(enable=not UseDegradationFactorInput));
    parameter Boolean UseDegradationFactorInput = false
      annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="External inputs/outputs"));
    protected
      Real DegFact;
    equation
      if not UseDegradationFactorInput then
        DegFact = DegradationFactorParam;
      end if;
    //  ProteinsMassConcentration2Concentration(c.u*1000) = q_in.conc;
      c.u = q_in.conc;
      q_in.q = DegradationBasic * c.val*DegFact;
    //  q_in.q =ProteinsMass2AmountOfSubstance(DegradationBasic*c.val*sec,ghostPlasmaVol)/sec;

    //TODO: state
    //der(degradedAmount) = q_in.q;
    //  ProteinsMass2AmountOfSubstance(degradedMass,ghostPlasmaVol) = degradedAmount;
     annotation (
        defaultComponentName="degradation",
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                100,100}}), graphics={Rectangle(
              extent={{-100,-50},{100,50}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid), Text(
              extent={{-88,-50},{100,50}},
              lineColor={0,0,255},
              textString="%name")}),  Diagram(coordinateSystem(preserveAspectRatio=true,
                       extent={{-100,-100},{100,100}}), graphics),
        Documentation(revisions="<html>
<p><i>2009-2010</i></p>
<p>Marek Matejak, Charles University, Prague, Czech Republic </p>
</html>"),        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics));
    end Degradation;

    model ProteinCharge

      Physiolibrary.Types.RealIO.pHInput pH
        annotation (Placement(transformation(extent={{-120,-110},{-80,-70}})));
      Physiolibrary.Chemical.Interfaces.ChemicalPort_a port_a
        annotation (Placement(transformation(extent={{80,-10},{100,10}})));
      Physiolibrary.Types.RealIO.ConcentrationInput AlbuminDifferenceMolarFlow
        annotation (Placement(transformation(extent={{-120,70},{-80,110}})));
      constant Real AlbMolarMass( final unit = "g/mol")= 66000;
      FullBloodAcidBase.FiggeFencl3 figgeFencl(
        pH = pH,
        pCO2=40,
        Pi=1.15,
        alb=AlbuminDifferenceMolarFlow*AlbMolarMass)
        annotation (Placement(transformation(extent={{-20,0},{0,20}})));
                //  parameter Real alb (unit = "g/dl")= 4.4;
      Real albuminChargeFlowRate( unit = "mol/s")= if UseSimpleAlbuminCharge then AlbuminDifferenceMolarFlow * AlbuminMoleCharge  else figgeFencl.atch*1000
        "In FF3 model, the units are in mEq/l concentration.";
      Real bicarbonateFlowRate( unit = "mol/s") = albuminChargeFlowRate
        "For each negative charge of albumin, a hydrogen ion is built, eats up one bicarbonate.";
       parameter Boolean UseSimpleAlbuminCharge = true;
        constant Real AlbuminMoleCharge = -18.6;
    equation
      port_a.q = -bicarbonateFlowRate "Ouflow, therefore negative";
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}})));
    end ProteinCharge;

    model AcidBaseBuffers

      Physiolibrary.Chemical.Interfaces.ChemicalPort_a port_a
        annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      Physiolibrary.Types.RealIO.pHOutput pH = figgeFencl3_1.pH
        annotation (Placement(transformation(extent={{-80,-100},{-100,-80}})));
      FiggeFencl3 figgeFencl3_1(
        SID=SID,
        pCO2=pCO2,
        Pi=1.45,
        alb=alb)
        annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
    /*  
  SAnomogram_formalization.SAVanSlyke sAVanSlyke
    annotation (Placement(transformation(extent={{0,20},{20,40}})));
*/
      Real BE;
      Real SID = NSID + BE;
      Real NSID;
      parameter Real pCO2 = 40;
      Real hco3(unit = "mmol/l") = figgeFencl3_1.HCO3*1000
        "mEq/l = mmol/l = mol/m3 = Eq/m3";
      Real hco3MM( unit="mol");
      Physiolibrary.Types.RealIO.ConcentrationInput albuminConcentration
        annotation (Placement(transformation(extent={{-120,70},{-80,110}})));

    parameter Boolean UseConstantAlb = true;
    constant Physiolibrary.Types.MolarMass AlbuminMolarMass = 66;
    Real alb = if UseConstantAlb then 4.4 else albuminConcentration*AlbuminMolarMass/10;
    parameter Physiolibrary.Types.Volume plasmaVol = 0.003;
    protected
      FiggeFencl3 normalPlasma(pH=7.4, pCO2 = 40, Pi = 1.15, alb = alb, SID = NSID)
        annotation (Placement(transformation(extent={{-58,0},{-38,20}})));

    equation
      // DEBUG
      when time > 1.5*24*3600 then
        reinit(hco3MM, 0);
      end when;

      port_a.conc = hco3;

    der(hco3MM) = port_a.q;
    BE = hco3MM/plasmaVol;

      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}})));
    end AcidBaseBuffers;

    model ProteinDivision "60% of total plasma protein mass are albumin"

      Physiolibrary.Types.RealIO.ConcentrationInput totalProteins
        annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
      Physiolibrary.Types.RealIO.ConcentrationOutput albumin
        annotation (Placement(transformation(extent={{90,50},{110,70}})));
      Physiolibrary.Types.RealIO.MassConcentrationOutput nonAlbumin
        annotation (Placement(transformation(extent={{90,-50},{110,-30}})));
      Physiolibrary.Types.RealIO.ConcentrationOutput nonAlbuminMolarConc
        annotation (Placement(transformation(extent={{90,-90},{110,-70}})));

    //  Physiolibrary.Types.MassConcentration tProtMassConc;
    //  Physiolibrary.Types.MassConcentration albuminMassConc;
    //  parameter Physiolibrary.Types.MolarMass albuminMM=66.5 "albumin molar mass";

    protected
      Physiolibrary.Types.MolarMass nonAlbuminMM
        "average molar mass of non-albumin proteins";

    equation
     /* totalProteins = ProteinsMassConcentration2Concentration(tProtMassConc);
  albuminMassConc = 0.6 * tProtMassConc;
  nonAlbumin = tProtMassConc - albuminMassConc;
  albumin = albuminMassConc / 66.5;
  nonAlbuminMolarConc = totalProteins - albumin; */

      albumin = totalProteins * (0.63/1.45);
      nonAlbuminMolarConc = totalProteins - albumin;
      nonAlbumin = nonAlbuminMM * nonAlbuminMolarConc;

    // inversion of totalProteins=(320*101325/760)/(310.15*8.314) *0.001*tProtMassConc + (1160*101325/760)/(310.15*8.314)* (0.001*tProtMassConc)^2;
    //  tProtMassConc = 0.0000170159 * (-8.106e6 + 63.6632 * ((1.6212e10 + 1.4208e10 * totalProteins)^0.5));
    // linear aproximation at point totalProteins = 1.45 mmol/l :
      nonAlbuminMM = 34.16-10*(totalProteins-1.45);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{100,100}})));
    end ProteinDivision;

    model pulse

      Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(
              extent={{80,-10},{100,10}}), iconTransformation(extent={{80,-10},{100,
                10}})));
      parameter Modelica.SIunits.Time startTime=0
        "Output = offset for time < startTime";
      parameter Modelica.SIunits.Time length=0
        "Output = offset for time < startTime";
      parameter Real dose;
      parameter Real normal;

    equation
      if time > startTime and time < startTime + length then
        y = dose;
      else
        y = normal;
      end if;

    end pulse;
  end AlbuminBorderFlux;

  package WBABB

    package Interfaces

      type BCE = enumeration(
          Water   "Water 1",
          K   "    K 2",
          Na   "   Na 3",
          Ca   "   Ca 4",
          Mg "Mg 5",
          Cl   "   Cl 6",
          Alb   "  Alb 7",
          Prot   " Other Proteins 8",
          P    "   P 9",
          tO2  "   tO2 10",
          tCO2 "   tCO2 11",
          BEox "   BEox 12",
          Hb   "   Hb 13 - concentration of Hemoglobin tetramer, i.e. consisting of four hems",
          Glucose "Glucose 14",
          Urea "Urea 15",
          Unchrg "Unchrg 16",
          Lactate "Lactate 17",
          Acetate "Acetate 18",
          Citrate "Citrate 19") "Blood contents enumeration";
      partial connector BloodConnector
        Real p( unit = "Pa");
        flow Real q( unit="m3/s");
        stream Real BC[BCE]( each unit="mmol/l")
                                                "Blood contents";
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Ellipse(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,0},
                fillColor={255,170,170},
                fillPattern=FillPattern.Solid)}),                      Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end BloodConnector;

      connector BloodIn
        extends Interfaces.BloodConnector;
        annotation (defaultComponentName = "bIn", Icon(graphics={Ellipse(extent={{-60,60},{60,-60}}, lineColor=
                    {0,0,0})}));
      end BloodIn;

      connector BloodOut
        extends Interfaces.BloodConnector;
        annotation (defaultComponentName = "bOut", Icon(graphics={Line(points={{
                    -60,60},{60,-60}},                                                             color={0,0,
                    0}),   Line(points={{60,60},{-60,-60}}, color={0,0,0})}));
      end BloodOut;
    end Interfaces;

    package Sources

      partial model BloodSourceBase
        import FullBloodAcidBase.WBABB.Interfaces.BCE;
        import FullBloodAcidBase.WBABB.Constants.*;

        Interfaces.BloodOut bloodOut
          annotation (Placement(transformation(extent={{80,80},{100,100}})));

      parameter Real K(unit = "mmol/l") = 4.4;
      parameter Real Water(unit = "mmol/l") = waterRatio*1000/waterMolarMass;
      parameter Real waterRatio(unit = "1") = 0.93;
      parameter Real Na(unit = "mmol/l")=140;
      parameter Real Ca(unit = "mmol/l")=1.5 "Ionized Ca";
      parameter Real Mg(unit = "mmol/l")=0.5 "Ionized Mg";
      parameter Real Cl(unit = "mmol/l")= 104;
      parameter Real Alb(unit = "mmol/l") = AlbMass / albuminMolarMass;
      parameter Real AlbMass( unit = "kg/m3") = 44;
      parameter Real Prot(unit = "mmol/l") = Alb/0.6*0.4;
      parameter Real P(unit = "mmol/l")=1.1;
      parameter Real tO2(unit = "mmol/l") = 8.8;
      parameter Real tCO2(unit = "mmol/l") = 21.7;
      parameter Real BEox(unit = "mmol/l")= 0;
      parameter Real Hb(unit = "mmol/l") = HbMass / HbMolarMass;
      parameter Real HbMass( unit = "kg/m3")=150;
      parameter Real Glucose(unit = "mmol/l") = 3;
      parameter Real Urea(unit = "mmol/l") = 3.6;
      parameter Real Unchrg(unit = "mmol/l") = 1;
      parameter Real Lactate(unit = "mmol/l") = 1;
      parameter Real Acetate(unit = "mmol/l") = 0;
      parameter Real Citrate(unit = "mmol/l") = 0;

      equation

        Water = bloodOut.BC[BCE.Water];
        K = bloodOut.BC[BCE.K];
        Na = bloodOut.BC[BCE.Na];
        Ca = bloodOut.BC[BCE.Ca];
        Mg = bloodOut.BC[BCE.Mg];
        Cl = bloodOut.BC[BCE.Cl];
        Alb = bloodOut.BC[BCE.Alb];
        Prot = bloodOut.BC[BCE.Prot];
        P = bloodOut.BC[BCE.P];
        tO2 = bloodOut.BC[BCE.tO2];
        tCO2 = bloodOut.BC[BCE.tCO2];
        BEox = bloodOut.BC[BCE.BEox];
        Hb = bloodOut.BC[BCE.Hb];
        Glucose = bloodOut.BC[BCE.Glucose];
        Urea = bloodOut.BC[BCE.Urea];
        Unchrg = bloodOut.BC[BCE.Unchrg];
        Lactate = bloodOut.BC[BCE.Lactate];
        Acetate = bloodOut.BC[BCE.Acetate];
        Citrate = bloodOut.BC[BCE.Citrate];

      /*  for i in BCE loop
    if i == BCE.Hb then
      bloodOut.BC[BCE.Hb] = 0;
    else
      bloodOut.BC[i] = 0;
    end if;
    end for;
*/
          /*
    Water   "1",
    K   "    2",
    Na   "   3",
    Ca   "   4",
    Cl   "   5",
    Alb   "  6",
    Prot   " 7",
    P    "   8",
    tO2  "   9",
    tCO2 "   10",
    BEox "   11",
    Hb   "   12",
    Unchrg " 13",
    Lactate,
    Acetate,
    Citrate) "Blood contents enumeration";
    
Hb = bloodOut.BC[BCE.Hb];
tO2 = bloodOut.BC[BCE.tO2];
tCO2 = bloodOut.BC[BCE.tCO2];
P = bloodOut.BC[BCE.P];
Alb = bloodOut.BC[BCE.Alb];
BEox = bloodOut.BC[BCE.BEox];
*/
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                                                     Text(
                extent={{-90,-110},{110,-10}},
                lineColor={0,0,127},
                fillColor={255,170,170},
                fillPattern=FillPattern.Solid,
                textString="%class"), Rectangle(extent={{-100,100},{100,-100}},
                  lineColor={28,108,200})}),                           Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end BloodSourceBase;

      model BloodFlowSource
        extends BloodSourceBase;
        parameter Real outFlow( unit = "m3/s", displayUnit = "ml/min");
      equation
        bloodOut.q + outFlow = 0;
      end BloodFlowSource;

      model BloodPressureSource
        extends BloodSourceBase;
        parameter Real pressure( unit = "Pa", displayUnit = "kPa");
      equation
        bloodOut.p = pressure;
      end BloodPressureSource;
    end Sources;

    package Constants
      final constant Real waterMolarMass( final unit = "kg/mol")= 0.01801528;
      final constant Real HbMolarMass( final unit = "kg/mol")= 64.500;
      final constant Real albuminMolarMass( final unit = "kg/mol")= 65.5;
      final constant Real pa2mmHg(final unit = "mmHg/Pa") = 1/133;
    end Constants;

    package Parts

      package Auxiliary

        partial model totalO2Base
          import FullBloodAcidBase.WBABB.Constants.*;
          Modelica.Blocks.Interfaces.RealInput pH
            annotation (Placement(transformation(extent={{60,-30},{100,10}})));
          Modelica.Blocks.Interfaces.RealInput Hb "concentration of hemoglobin tetramer"
            annotation (Placement(transformation(extent={{-120,-110},{-80,-70}})));
          Modelica.Blocks.Interfaces.RealInput tO2
            annotation (Placement(transformation(extent={{-120,10},{-80,50}})));
          Modelica.Blocks.Interfaces.RealInput pCO2
            annotation (Placement(transformation(extent={{-20,-110},{20,-70}})));

          Modelica.Blocks.Interfaces.RealOutput pO2
            annotation (Placement(transformation(extent={{80,80},{100,100}})));
          Modelica.Blocks.Interfaces.RealOutput sO2
            annotation (Placement(transformation(extent={{40,-100},{60,-80}})));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));
        /*
Real Hbgdl(unit = "g/dl") = Hb * HbMolarMass/10 "Hemoglobin in g per deciliter";
Real pCO2mmHg( unit = "mmHg") = pCO2*pa2mmHg "pCO2 in torr";
Real pO2mmHg( unit = "mmHg") = pO2*pa2mmHg "pO2 in torr";
*/
        end totalO2Base;

        partial model totalCO2Base
          import FullBloodAcidBase.WBABB.Constants.*;

          Modelica.Blocks.Interfaces.RealInput pH
            annotation (Placement(transformation(extent={{60,70},{100,110}})));
          Modelica.Blocks.Interfaces.RealInput Hb
            annotation (Placement(transformation(extent={{-120,-30},{-80,10}})));
          Modelica.Blocks.Interfaces.RealInput tCO2
            annotation (Placement(transformation(extent={{-120,-70},{-80,-30}})));
          Modelica.Blocks.Interfaces.RealInput sO2
            annotation (Placement(transformation(extent={{60,10},{100,50}})));
          Modelica.Blocks.Interfaces.RealInput pO2
            annotation (Placement(transformation(extent={{60,-50},{100,-10}})));
          Modelica.Blocks.Interfaces.RealOutput pCO2( start=6000)
            annotation (Placement(transformation(extent={{-100,70},{-60,110}})));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));
        end totalCO2Base;

        partial model ABBBase

          Modelica.Blocks.Interfaces.RealInput Hb
            annotation (Placement(transformation(extent={{-120,70},{-80,110}})));
          Modelica.Blocks.Interfaces.RealInput P
            annotation (Placement(transformation(extent={{-120,-30},{-80,10}})));
          Modelica.Blocks.Interfaces.RealInput Alb
            annotation (Placement(transformation(extent={{-120,-70},{-80,-30}})));
          Modelica.Blocks.Interfaces.RealInput BEox
            annotation (Placement(transformation(extent={{-120,-110},{-80,-70}})));
          Modelica.Blocks.Interfaces.RealInput pCO2
            annotation (Placement(transformation(extent={{-120,80},{-80,40}})));
          Modelica.Blocks.Interfaces.RealOutput pH
            annotation (Placement(transformation(extent={{80,70},{120,110}})));
          Modelica.Blocks.Interfaces.RealOutput HCO3
            annotation (Placement(transformation(extent={{80,30},{120,70}})));
          Modelica.Blocks.Interfaces.RealInput sO2
            annotation (Placement(transformation(extent={{-120,50},{-80,10}})));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));
        end ABBBase;

        model ABBVanSlyke
          extends ABBBase;
          import FullBloodAcidBase.WBABB.Constants.*;
          SAnomogram_formalization.SAVanSlyke sAVanSlyke(
            pCO2=pCO2/133,
            BEox=BEox,
            Hct=Hb  * HbMolarMass/10/ 33.34,
            sO2=sO2)
            annotation (Placement(transformation(extent={{-20,0},{0,20}})));
        equation
          HCO3 = sAVanSlyke.cHCO3;
          pH = sAVanSlyke.pH;
        end ABBVanSlyke;

        model totalO2Physiomodel
          extends totalO2Base;

          Real aO2 = exp(log(0.0105)+(-0.0115*(T-T0))+0.5*0.00042*(T-T0)^2)/1000 "o2 solubility";
          Physiolibrary.Types.Fraction sO2CO(start=0.75) "What is this?";
          Physiolibrary.Types.Pressure pO2CO;
          Physiolibrary.Types.Concentration cO2Hb(start=6);
          //
          Physiolibrary.Types.Fraction sCO;

          Physiolibrary.Types.Concentration ceHb = ctHb * (1-FCOHb-FMetHb) "effective haemoglobin";
          Real a(start=0.5) = dadpH*(pH-pH0)+dadlnpCO2*log(max(1e-15+1e-22*pCO2,pCO2/pCO20)) +dadxMetHb*FMetHb+(dadcDPG0 + dadcDPGxHbF*FHbF)*(cDPG/cDPG0 - 1);

          Real x=log(pO2CO/7000) - a - 0.055*(T-T0); //namiesto:  x=log(pO2CO/7) - a - 0.055*(T-37);;
          Real y;
          Real k=0.5342857;
          Real h = 3.5 + a;

         Physiolibrary.Types.Fraction FCOHb(start=0);
         Real ctHb = Hb*4 "Concentration of hemoglobin monomer";
         parameter Real T(start=310.15) = 273 + 37;
         parameter Real pCO = 0;
         parameter Real cDPG = 5;

         parameter Real FMetHb = 0;

         parameter Real FHbF = 0;

         Physiolibrary.Types.Concentration cdO2 = aO2*pO2;

         parameter Physiolibrary.Types.Temperature T0 = 273.15+37
            "normal temperature";
         parameter Physiolibrary.Types.pH pH0 = 7.4 "normal arterial pH";
         parameter Physiolibrary.Types.Pressure pCO20 = 5330
            "normal arterial CO2 partial pressure";
         parameter Real cDPG0 = 5
            "normal DPG,used by a";
         parameter Real dadcDPG0 = 0.3 "used by a";
         parameter Real dadcDPGxHbF = -0.1 "or perhabs -0.125";
         parameter Real dadpH = -0.88 "used by a";
         parameter Real dadlnpCO2 = 0.048 "used by a";
         parameter Real dadxMetHb = -0.7 "used by a";
         parameter Real dadxHbF = -0.25 "used by a";

        equation

          assert(tO2 <= ceHb*(1.06), "Model does not support this high level of oxygen in blood. Maximum of oxygen concentration should be connected with efective hemoglobin concentration!");
          tO2 = aO2*pO2 + ceHb*sO2;
          sO2 = cO2Hb/ceHb;

          //orginal:
          y-1.8747=x+h*tanh(k*x);
          y=log(sO2CO/(1-sO2CO));

            {pCO,FCOHb,pO2CO,sO2CO}=homotopy({sCO*pO2CO/ 218*sO2CO,sCO*(1-FMetHb),pO2 + 218*pCO,(cO2Hb + ctHb*FCOHb)/(ctHb*(1-FMetHb))},
            {0,0,pO2,sO2});

        end totalO2Physiomodel;

        model totalCO2Physiomodel
          extends totalCO2Base;
          import Modelica.Math;
          final constant Real T = 273 + 37;
          Physiolibrary.Types.Concentration tCO2_P(start=24, displayUnit="mmol/l") = cHCO3 + cdCO2;
          Real pK_ery = 6.125 - log10(1+10^(pH_ery-7.84-0.06*sO2));
          Physiolibrary.Types.GasSolubility aCO2_ery( displayUnit="mmol/l/mmHg") = 0.000195; //solubility 0.23 (mmol/l)/kPa at 25degC
          Physiolibrary.Types.Concentration tCO2_ery( displayUnit="mmol/l")= aCO2_ery*pCO2*(1+10^(pH_ery-pK_ery));
          Real ctHb = Hb*4 "Concentration of hemoglobin monomer";
          Real pH_ery = homotopy(7.19 + 0.77*(pH-7.4) + 0.035*(1-sO2),7.19 + 0.77*(pH-7.4))  "outgoing intracellular erytrocytes pH";
          Real Hct = (ctHb*0.44)/8.4 "hematocrit (erytrocytes volume/blood volume)";
          Real pK = 6.1 + (-0.0026)*(T-310.15);
          Real aCO2(final displayUnit="mmol/(l.kPa)") = 0.00023 * 10^(-0.0092*(T-310.15)); //solubility depends on temperature;
          Physiolibrary.Types.Concentration cdCO2(displayUnit="mmol/l") = aCO2*pCO2;
          Real cHCO3 = cdCO2 * 10^(pH-pK);

        equation

          tCO2 = tCO2_ery*Hct + tCO2_P*(1-Hct);

        end totalCO2Physiomodel;

        partial model OneBloodPort

          Interfaces.BloodIn bIn "Blood Inflow connector"
            annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
          Interfaces.BloodOut bOut "Blood Outflow connector"
            annotation (Placement(transformation(extent={{80,-12},{100,8}})));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                  Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200})}),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end OneBloodPort;
      end Auxiliary;

      model Tissues
        extends Auxiliary.OneBloodPort;
        import FullBloodAcidBase.WBABB.Interfaces.BCE;
        parameter Real resistance;
        parameter Real consumedO2(unit="mol/s", displayUnit = "mmol/s");
        parameter Real producedCO2(unit="mol/s", displayUnit = "mmol/s");
        //  parameter Real volume;
      equation
        bIn.q + bOut.q = 0;
        bIn.p = bOut.q*resistance;

        for i in BCE loop
          if i == BCE.tO2 then
            bOut.BC[BCE.tO2] = inStream(bIn.BC[BCE.tO2]) - consumedO2/bIn.q;
          elseif i == BCE.tCO2 then
            bOut.BC[BCE.tCO2] = inStream(bIn.BC[BCE.tCO2]) + producedCO2/bIn.q;
          else
            bOut.BC[i] = inStream(bIn.BC[i]);
          end if;
        end for;


        assert(bIn.q >= 0, "The blood in tissues cannot flow backwards!");
        for i in BCE loop
            0 = bIn.BC[i] "This should never happen, as the flow is projected from In to Out only";
        end for;

        annotation (Icon(graphics={Line(
                points={{-76,0},{72,0}},
                color={28,108,200},
                arrow={Arrow.None,Arrow.Filled}),                    Text(
                extent={{-100,-120},{100,-20}},
                lineColor={0,0,127},
                fillColor={255,170,170},
                fillPattern=FillPattern.Solid,
                textString="%class")}));
      end Tissues;

      model bloodABBMeasurement
        import FullBloodAcidBase.WBABB.Interfaces.BCE;

        Interfaces.BloodIn bloodIn
          annotation (Placement(transformation(extent={{80,80},{100,100}}),
              iconTransformation(extent={{80,80},{100,100}})));

        replaceable Auxiliary.totalO2Physiomodel
                                        totalO2Base1 constrainedby
          Auxiliary.totalO2Base
          annotation (Placement(transformation(extent={{-30,48},{-10,68}})));
        replaceable Auxiliary.totalCO2Physiomodel totalCO2Base1 constrainedby
          Auxiliary.totalCO2Base
          annotation (Placement(transformation(extent={{-20,-50},{0,-30}})));
        replaceable Auxiliary.ABBVanSlyke aBBB1 constrainedby Auxiliary.ABBBase
          annotation (Placement(transformation(extent={{18,-8},{50,24}})));
      protected
        Modelica.Blocks.Interfaces.RealInput Hb = inStream(bloodIn.BC[BCE.Hb])
          annotation (Placement(transformation(extent={{-120,2},{-80,42}})));
        Modelica.Blocks.Interfaces.RealInput tO2 = inStream(bloodIn.BC[BCE.tO2])
          annotation (Placement(transformation(extent={{-120,40},{-80,80}})));
        Modelica.Blocks.Interfaces.RealInput tCO2 = inStream(bloodIn.BC[BCE.tCO2])
          annotation (Placement(transformation(extent={{-120,-66},{-80,-26}})));
        Modelica.Blocks.Interfaces.RealInput P = inStream(bloodIn.BC[BCE.P])
          annotation (Placement(transformation(extent={{-8,2},{2,12}})));
        Modelica.Blocks.Interfaces.RealInput Alb = inStream(bloodIn.BC[BCE.Alb])
          annotation (Placement(transformation(extent={{-8,-6},{2,4}})));
        Modelica.Blocks.Interfaces.RealInput BEox = inStream(bloodIn.BC[BCE.BEox])
          annotation (Placement(transformation(extent={{-8,-14},{2,-4}})));
      equation
       bloodIn.q = 0;
       bloodIn.BC = zeros(size(bloodIn.BC, 1)) "nothing goes out";

        connect(aBBB1.pH, totalO2Base1.pH) annotation (Line(points={{50,22.4},{48,22.4},
                {48,22},{54,22},{54,58},{-14,58},{-14,57},{-12,57}},
                                               color={0,0,127}));
        connect(totalO2Base1.Hb, Hb) annotation (Line(points={{-30,49},{-32,49},{-32,22},
                {-100,22}}, color={0,0,127}));
        connect(aBBB1.Hb, Hb) annotation (Line(points={{18,22.4},{16,22.4},{16,22},{-100,
                22}},      color={0,0,127}));
        connect(totalCO2Base1.Hb, Hb) annotation (Line(points={{-20,-41},{-32,-41},{-32,
                22},{-100,22}}, color={0,0,127}));
        connect(aBBB1.pH, totalCO2Base1.pH) annotation (Line(points={{50,22.4},{48,22.4},
                {48,22},{54,22},{54,-30},{-2,-30},{-2,-31}},
                                            color={0,0,127}));
        connect(totalO2Base1.tO2, tO2) annotation (Line(points={{-30,61},{-100,61},{-100,
                62},{-100,62},{-100,60},{-100,60}},
                               color={0,0,127}));
        connect(totalO2Base1.pCO2, totalCO2Base1.pCO2) annotation (Line(points={{-20,49},
                {-18,49},{-18,20},{-18,-31}},          color={0,0,127}));
        connect(totalO2Base1.sO2, totalCO2Base1.sO2) annotation (Line(points={{-15,49},
                {-8,49},{-8,-36},{-8,-37},{-2,-37}}, color={0,0,127}));
        connect(totalCO2Base1.pO2, totalO2Base1.pO2) annotation (Line(points={{-2,-43},
                {-2,-43},{-2,-44},{62,-44},{62,67},{-11,67}}, color={0,0,127}));
        connect(totalCO2Base1.pCO2,aBBB1. pCO2) annotation (Line(points={{-18,-31},{-18,
                4},{-18,17.6},{18,17.6}},
                                    color={0,0,127}));
        connect(totalO2Base1.sO2, aBBB1.sO2)
          annotation (Line(points={{-15,49},{-8,49},{-8,14},{18,14},{18,12.8}},
                                                              color={0,0,127}));
        connect(aBBB1.P, P) annotation (Line(points={{18,6.4},{18,6},{18,6},{16,6},{16,
                7},{-3,7}},
              color={0,0,127}));
        connect(aBBB1.Alb, Alb)
          annotation (Line(points={{18,0},{18,0},{18,2},{18,2},{18,-1},{-3,-1}},
                                                                   color={0,0,127}));
        connect(BEox, aBBB1.BEox) annotation (Line(points={{-3,-9},{18,-9},{18,-6},{18,
                -6},{18,-6},{18,-6.4}},
                            color={0,0,127}));
        connect(totalCO2Base1.tCO2, tCO2) annotation (Line(points={{-20,-45},{-20,-45},
                {-20,-46},{-100,-46}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics
              ={                                                     Text(
                extent={{-100,-120},{100,-20}},
                lineColor={0,0,127},
                fillColor={255,170,170},
                fillPattern=FillPattern.Solid,
                textString="%class"), Rectangle(extent={{-100,100},{100,-100}},
                  lineColor={0,255,0})}),                              Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end bloodABBMeasurement;
    end Parts;

    package Tests

      model TissuesTest

        Sources.BloodFlowSource bloodFlowSource(
          HbMass=150,
          Hb=9.3,
          outFlow=8.3333333333333e-05)
          annotation (Placement(transformation(extent={{-82,12},{-64,32}})));
        Sources.BloodPressureSource bloodPressureSource(pressure=0)
          annotation (Placement(transformation(extent={{20,12},{40,32}})));
        Parts.Tissues tissues(
          resistance=1,
          consumedO2(displayUnit="mmol/day") = 0.00025462962962963,
          producedCO2(displayUnit="mmol/day") = 0.00025462962962963)
          annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
        Parts.bloodABBMeasurement bloodABBMeasurementTest
          annotation (Placement(transformation(extent={{-74,-20},{-54,0}})));
        Parts.bloodABBMeasurement bloodABBMeasurementTest1
          annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
      equation
        connect(bloodFlowSource.bloodOut, tissues.bIn) annotation (Line(points={{-64.9,
                31},{-64.9,30},{-39,30}},     color={0,0,0}));
        connect(bloodPressureSource.bloodOut, tissues.bOut) annotation (Line(
              points={{39,31},{4,31},{4,29.8},{-21,29.8}}, color={0,0,0}));
        connect(tissues.bIn, bloodABBMeasurementTest.bloodIn) annotation (Line(
              points={{-39,30},{-44,30},{-44,-1},{-55,-1}}, color={0,0,0}));
        connect(tissues.bOut, bloodABBMeasurementTest1.bloodIn) annotation (
            Line(points={{-21,29.8},{-12,29.8},{-12,-1},{-1,-1}}, color={0,0,0}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TissuesTest;

      model ABRtest
        Parts.Auxiliary.ABBVanSlyke aBBVanSlyke
          annotation (Placement(transformation(extent={{-10,-64},{92,38}})));
        Modelica.Blocks.Sources.Constant const(k=2.32)
          annotation (Placement(transformation(extent={{-98,32},{-78,52}})));
        Modelica.Blocks.Sources.Constant const1(k=5400)
          annotation (Placement(transformation(extent={{-96,-4},{-76,16}})));
        Modelica.Blocks.Sources.Constant const2(k=1)
          annotation (Placement(transformation(extent={{-100,-38},{-80,-18}})));
        Modelica.Blocks.Sources.Constant const3(k=1.1)
          annotation (Placement(transformation(extent={{-98,-72},{-78,-52}})));
        Modelica.Blocks.Sources.Constant const4(k=0.67) annotation (Placement(
              transformation(extent={{-100,-102},{-80,-82}})));
        Modelica.Blocks.Sources.Constant const5(k=0)
          annotation (Placement(transformation(extent={{-32,-100},{-12,-80}})));
      equation
        connect(aBBVanSlyke.Hb, const.y) annotation (Line(points={{-10,32.9},{
                -40,32.9},{-40,42},{-77,42}}, color={0,0,127}));
        connect(aBBVanSlyke.pCO2, const1.y) annotation (Line(points={{-10,17.6},
                {-43,17.6},{-43,6},{-75,6}}, color={0,0,127}));
        connect(aBBVanSlyke.sO2, const2.y) annotation (Line(points={{-10,2.3},{
                -44,2.3},{-44,-28},{-79,-28}}, color={0,0,127}));
        connect(aBBVanSlyke.P, const3.y) annotation (Line(points={{-10,-18.1},{
                -43,-18.1},{-43,-62},{-77,-62}}, color={0,0,127}));
        connect(aBBVanSlyke.Alb, const4.y) annotation (Line(points={{-10,-38.5},
                {-42,-38.5},{-42,-92},{-79,-92}}, color={0,0,127}));
        connect(aBBVanSlyke.BEox, const5.y) annotation (Line(points={{-10,-58.9},
                {-10,-73.45},{-11,-73.45},{-11,-90}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end ABRtest;

      model CO2Test
        Parts.Auxiliary.totalCO2Physiomodel totalCO2SA
          annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
        Modelica.Blocks.Sources.Constant const(k=9)
          annotation (Placement(transformation(extent={{-80,-18},{-60,2}})));
        Modelica.Blocks.Sources.Constant const1(k=25)
          annotation (Placement(transformation(extent={{-82,-52},{-62,-32}})));
        Modelica.Blocks.Sources.Constant const2(k=7.4)
          annotation (Placement(transformation(extent={{-40,64},{-20,84}})));
        Modelica.Blocks.Sources.Constant const3(k=1)
          annotation (Placement(transformation(extent={{36,-6},{56,14}})));
        Modelica.Blocks.Sources.Constant const4(k=10)
          annotation (Placement(transformation(extent={{30,-52},{50,-32}})));
      equation
        connect(totalCO2SA.Hb, const.y) annotation (Line(points={{-20,-11},{-40,
                -11},{-40,-8},{-59,-8}}, color={0,0,127}));
        connect(totalCO2SA.tCO2, const1.y) annotation (Line(points={{-20,-15},{
                -41,-15},{-41,-42},{-61,-42}}, color={0,0,127}));
        connect(const2.y, totalCO2SA.pH) annotation (Line(points={{-19,74},{-10,
                74},{-10,-1},{-2,-1}}, color={0,0,127}));
        connect(const3.y, totalCO2SA.sO2) annotation (Line(points={{57,4},{28,4},
                {28,-7},{-2,-7}}, color={0,0,127}));
        connect(const4.y, totalCO2SA.pO2) annotation (Line(points={{51,-42},{24,
                -42},{24,-13},{-2,-13}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end CO2Test;

      model O2Test
        Parts.Auxiliary.totalO2Physiomodel totalO2SA
          annotation (Placement(transformation(extent={{-20,-18},{0,2}})));
        Modelica.Blocks.Sources.Constant const(k=8)
          annotation (Placement(transformation(extent={{-86,-4},{-66,16}})));
        Modelica.Blocks.Sources.Constant const1(k=8)
          annotation (Placement(transformation(extent={{-78,-36},{-58,-16}})));
        Modelica.Blocks.Sources.Constant const2(k=5300)
          annotation (Placement(transformation(extent={{-34,-64},{-14,-44}})));
        Modelica.Blocks.Sources.Constant const3(k=7.4)
          annotation (Placement(transformation(extent={{24,-26},{44,-6}})));
      equation
        connect(const.y, totalO2SA.tO2) annotation (Line(points={{-65,6},{-44,6},
                {-44,-5},{-20,-5}}, color={0,0,127}));
        connect(totalO2SA.Hb, const1.y) annotation (Line(points={{-20,-17},{-36,
                -17},{-36,-26},{-57,-26}}, color={0,0,127}));
        connect(totalO2SA.pH, const3.y) annotation (Line(points={{-2,-9},{22,-9},
                {22,-16},{45,-16}}, color={0,0,127}));
        connect(totalO2SA.pCO2, const2.y) annotation (Line(points={{-10,-17},{
                -13,-17},{-13,-54}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end O2Test;
    end Tests;

    package Thrash

      model totalO2Physiomodel2
        extends Parts.Auxiliary.totalO2Base;

        Real aO2 = exp(log(0.0105)+(-0.0115*(T-T0))+0.5*0.00042*(T-T0)^2)/1000 "o2 solubility";
        Physiolibrary.Types.Concentration cO2Hb(start=6) = ceHb*sO2;
        Physiolibrary.Types.Concentration ceHb = ctHb "effective haemoglobin";
        Real a(start=0.5) = dadpH*(pH-pH0)+dadlnpCO2*log(max(1e-15+1e-22*pCO2,pCO2/pCO20)) +(dadcDPG0)*(cDPG/cDPG0 - 1);
        Real x=log(pO2/7000) - a - 0.055*(T-T0); //namiesto:  x=log(pO2CO/7) - a - 0.055*(T-37);;
        Real y= x+h*tanh(k*x) + 1.8747;
        Real k=0.5342857;
        Real h = 3.5 + a;

       Real ctHb = Hb;

       Physiolibrary.Types.Concentration cdO2 = aO2*pO2;

       parameter Real T(start=310.15) = 273 + 37;
       parameter Real cDPG = 5;
       parameter Physiolibrary.Types.Temperature T0 = 273.15+37
          "normal temperature";
       parameter Physiolibrary.Types.pH pH0 = 7.4 "normal arterial pH";
       parameter Physiolibrary.Types.Pressure pCO20 = 5330
          "normal arterial CO2 partial pressure";
       parameter Real cDPG0 = 5
          "normal DPG,used by a";
       parameter Real dadcDPG0 = 0.3 "used by a";
       parameter Real dadcDPGxHbF = -0.1 "or perhabs -0.125";
       parameter Real dadpH = -0.88 "used by a";
       parameter Real dadlnpCO2 = 0.048 "used by a";
       parameter Real dadxMetHb = -0.7 "used by a";
       parameter Real dadxHbF = -0.25 "used by a";
       Real eh = (cO2Hb)/(ctHb);
      equation

        assert(tO2 <= ceHb*(1.06), "Model does not support this high level of oxygen in blood. Maximum of oxygen concentration should be connected with efective hemoglobin concentration!");
        tO2 = aO2*pO2 + ceHb*sO2;
        sO2 = log(eh/(1-eh));
       // sO2 = cO2Hb/ceHb;

      end totalO2Physiomodel2;

      model totalCO2SA
        extends Parts.Auxiliary.totalCO2Base;
        final constant Real T(final unit = "K") = 273.15 + 37;

        Real pK=6.1 + (-0.0026)*(T - 310.15)
          "Henderson-Hasselbalch equation: pK for HCO3";
        Real aCO2=0.00023*10^(-0.0092*(T - 310.15))
          "solubility of CO2 depends on temperature";

        //CO2
        Real cFreeCO2=aCO2*pCO2 "Dissolved CO2";
        Real cHCO3=cFreeCO2*10^(pH - pK) "dissolved HCO3";
      equation
          tCO2 = cHCO3 + cFreeCO2;
      end totalCO2SA;

      model totalO2SA
        "Calculations of total O2 and pO2 according to Siggaard-Andersen 54 variables.."
        extends Parts.Auxiliary.totalO2Base;

       /*  
   // O2
  // interchangeable pO2 / tO2 input or output
  constant Real FMetHb=0;
  constant Real FCOHb=0;
  constant Real cDPG=5;
  final constant Real T(final unit = "K") = 273.15 + 37;

//   Real a1=-0.72*(pH - 7.4);
//   Real a2=0.09*log(pCO2/5330);
//   Real a3=0.7*FMetHb;
//   Real a4=(0.3 - (0.1*0.005))*((cDPG/5) - 1);
//   Real a5=-0.25*0.005;
//  Real a=a1 + a2 + a3 + a4 + a5;
  Real a=-0.88*(pH-7.4)+0.048*log(max(1e-15+1e-19*pCO2,pCO2/5.33))-0.7*FMetHb+(0.3-0.25*FCOHb)*cDPG/(5-1); //Bohr coefficient: -der(log10(pO2),pH)=0.88/ln(10)=0.38, -der(ln(pO2),pH)=0.88
//  Real k0=0.5343;
//  Real b=0.055*(T - 37 - 273);
 Real y=log(0.867/(1 - 0.867));
//  Real h1=3.5 + a;
 // Real x0=a + b;
  Real cHb=(1 - (FMetHb + FCOHb))*Hb;
  Real O2Solubility = exp(log(0.0105)+(-0.0115*(T-37-273))+0.5*0.00042*(T-37-273)^2)/1000; //solubility

  //exp(log(0.0000105) + (-0.0000115*(T - 37 - 273)) + 2*0.000000105    *(T - 37 - 273)^2);
  // http://www.siggaard-andersen.dk/OsaTextbook.htm
  Real p=pO2 + (218*0.00023);
//  Real x1=log(pO2 + (218*0.00023)/7000);
  //Real y=y0 + (x1 - x0) + (h1*tanh(k0*(x1 - x0)));

  Real x=log(pO2/7000) - a - 0.055*(T-37-273); //namiesto:  x=log(pO2CO/7) - a - 0.055*(T-37);
//  Real y =x+h*tanh(k*x) + 1.8747;
  Real k=0.5342857;
  Real h=3.5 + a;
//  y=log(sO2CO/(1-sO2CO));

  Real s=exp(y)/(1 + exp(y));
  Real sO2_=(Hb*(s*(1 - FMetHb) - FCOHb))/cHb;
  Real freeO2=O2Solubility*pO2;
  Real tO2_=freeO2 + sO2*(Hb - (FMetHb*Hb) - (FCOHb*Hb));
  Real Hct=(cHb*0.44)/8.4;
  */

              constant Real FMetHb=0;
              constant Real FCOHb=0;
              constant Real T = 273.15 + 37;
              constant Real cDPG=5 "TODO";

              Real a1=-0.88*(pH - 7.4);
              Real a2=0.048*log(pCO2/5300);
              Real a3=-0.7*FMetHb;
              Real a4=(0.3 - (0.1*0.005))*((cDPG/5) - 1);
              Real a5=-0.25*0.005;
              Real a=a1 + a2 + a3 + a4 + a5;
              Real k0=0.5343;
              Real b=0.055*(T - 37 - 273);
              Real y0=log(0.867/(1 - 0.867));
              Real h1=3.5 + a;
              Real x0=a + b;
              Real cHb=(1 - (FMetHb + FCOHb))*Hb;
              Real O2Solubility=exp(log(0.0000105) + (-0.0000115*(T - 37 - 273)) + 2*0.000000105
                  *(T - 37 - 273)^2);
              // http://www.siggaard-andersen.dk/OsaTextbook.htm
              Real p=pO2 + (218*0.00023);
              Real x1=log(pO2 + (218*0.00023)/7000);
              Real y=y0 + (x1 - x0) + (h1*tanh(k0*(x1 - x0)));
              Real s=exp(y)/(1 + exp(y));
              Real sO2_=(Hb*(s*(1 - FMetHb) - FCOHb))/cHb;
              Real freeO2=O2Solubility*pO2;
              Real tO2_=freeO2 + sO2*(Hb - (FMetHb*Hb) - (FCOHb*Hb));
              Real Hct=(cHb*0.44)/8.4;
      equation
        //y =x+h*tanh(k*x) + 1.8747;
        sO2_ = sO2;
        tO2_ = tO2;
      end totalO2SA;

      model bloodABBMeasurementTest
        import FullBloodAcidBase.WBABB.Interfaces.BCE;

        Interfaces.BloodIn bloodIn
          annotation (Placement(transformation(extent={{80,80},{100,100}}),
              iconTransformation(extent={{80,80},{100,100}})));

        replaceable Thrash.totalCO2SA totalCO2Base1 constrainedby
          Parts.Auxiliary.totalCO2Base
          annotation (Placement(transformation(extent={{-20,-50},{0,-30}})));
        replaceable Parts.Auxiliary.ABBVanSlyke aBBB1 constrainedby
          Parts.Auxiliary.ABBBase
          annotation (Placement(transformation(extent={{18,-8},{50,24}})));
      protected
        Modelica.Blocks.Interfaces.RealInput Hb = inStream(bloodIn.BC[BCE.Hb])
          annotation (Placement(transformation(extent={{-120,2},{-80,42}})));
        Modelica.Blocks.Interfaces.RealInput tO2 = inStream(bloodIn.BC[BCE.tO2])
          annotation (Placement(transformation(extent={{-120,40},{-80,80}})));
        Modelica.Blocks.Interfaces.RealInput tCO2 = inStream(bloodIn.BC[BCE.tCO2])
          annotation (Placement(transformation(extent={{-120,-66},{-80,-26}})));
        Modelica.Blocks.Interfaces.RealInput P = inStream(bloodIn.BC[BCE.P])
          annotation (Placement(transformation(extent={{-8,2},{2,12}})));
        Modelica.Blocks.Interfaces.RealInput Alb = inStream(bloodIn.BC[BCE.Alb])
          annotation (Placement(transformation(extent={{-8,-6},{2,4}})));
        Modelica.Blocks.Interfaces.RealInput BEox = inStream(bloodIn.BC[BCE.BEox])
          annotation (Placement(transformation(extent={{-8,-14},{2,-4}})));
      public
        Modelica.Blocks.Sources.Constant const(k=1)
          annotation (Placement(transformation(extent={{-40,64},{-20,84}})));
        Modelica.Blocks.Sources.Constant const1(k=10)
          annotation (Placement(transformation(extent={{-30,-92},{-10,-72}})));
      equation
       bloodIn.q = 0;
       bloodIn.BC = zeros(size(bloodIn.BC, 1)) "nothing goes out";

        connect(aBBB1.Hb, Hb) annotation (Line(points={{18,22.4},{16,22.4},{16,22},{-100,
                22}},      color={0,0,127}));
        connect(totalCO2Base1.Hb, Hb) annotation (Line(points={{-20,-41},{-32,-41},{-32,
                22},{-100,22}}, color={0,0,127}));
        connect(aBBB1.pH, totalCO2Base1.pH) annotation (Line(points={{50,22.4},{48,22.4},
                {48,22},{54,22},{54,-30},{-2,-30},{-2,-31}},
                                            color={0,0,127}));
        connect(totalCO2Base1.pCO2,aBBB1. pCO2) annotation (Line(points={{-18,-31},{-18,
                4},{-18,17.6},{18,17.6}},
                                    color={0,0,127}));
        connect(aBBB1.P, P) annotation (Line(points={{18,6.4},{18,6},{18,6},{16,6},{16,
                7},{-3,7}},
              color={0,0,127}));
        connect(aBBB1.Alb, Alb)
          annotation (Line(points={{18,0},{18,0},{18,2},{18,2},{18,-1},{-3,-1}},
                                                                   color={0,0,127}));
        connect(BEox, aBBB1.BEox) annotation (Line(points={{-3,-9},{18,-9},{18,-6},{18,
                -6},{18,-6},{18,-6.4}},
                            color={0,0,127}));
        connect(totalCO2Base1.tCO2, tCO2) annotation (Line(points={{-20,-45},{-20,-45},
                {-20,-46},{-100,-46}}, color={0,0,127}));
        connect(const.y, aBBB1.sO2) annotation (Line(points={{-19,74},{8,74},{8,
                12.8},{18,12.8}}, color={0,0,127}));
        connect(const.y, totalCO2Base1.sO2) annotation (Line(points={{-19,74},{
                -10,74},{-10,-22},{10,-22},{10,-37},{-2,-37}}, color={0,0,127}));
        connect(const1.y, totalCO2Base1.pO2) annotation (Line(points={{-9,-82},
                {0,-82},{0,-78},{16,-78},{16,-43},{-2,-43}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics
              ={                                                     Text(
                extent={{-100,-120},{100,-20}},
                lineColor={0,0,127},
                fillColor={255,170,170},
                fillPattern=FillPattern.Solid,
                textString="%class"), Rectangle(extent={{-100,100},{100,-100}},
                  lineColor={0,255,0})}),                              Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end bloodABBMeasurementTest;

      model bloodABBMeasurementTestHack
        import FullBloodAcidBase.WBABB.Interfaces.BCE;

        Interfaces.BloodIn bloodIn
          annotation (Placement(transformation(extent={{80,80},{100,100}}),
              iconTransformation(extent={{80,80},{100,100}})));

        replaceable Parts.Auxiliary.totalO2Physiomodel totalO2Base1
          constrainedby Parts.Auxiliary.totalO2Base
          annotation (Placement(transformation(extent={{-30,48},{-10,68}})));
        replaceable Parts.Auxiliary.totalCO2Physiomodel totalCO2Base1
          constrainedby Parts.Auxiliary.totalCO2Base
          annotation (Placement(transformation(extent={{-20,-50},{0,-30}})));
        replaceable Parts.Auxiliary.ABBVanSlyke aBBB1 constrainedby
          Parts.Auxiliary.ABBBase
          annotation (Placement(transformation(extent={{18,-8},{50,24}})));
      protected
        Modelica.Blocks.Interfaces.RealInput Hb = inStream(bloodIn.BC[BCE.Hb])
          annotation (Placement(transformation(extent={{-120,2},{-80,42}})));
        Modelica.Blocks.Interfaces.RealInput tO2 = inStream(bloodIn.BC[BCE.tO2])
          annotation (Placement(transformation(extent={{-120,40},{-80,80}})));
        Modelica.Blocks.Interfaces.RealInput tCO2 = inStream(bloodIn.BC[BCE.tCO2])
          annotation (Placement(transformation(extent={{-120,-66},{-80,-26}})));
        Modelica.Blocks.Interfaces.RealInput P = inStream(bloodIn.BC[BCE.P])
          annotation (Placement(transformation(extent={{-8,2},{2,12}})));
        Modelica.Blocks.Interfaces.RealInput Alb = inStream(bloodIn.BC[BCE.Alb])
          annotation (Placement(transformation(extent={{-8,-6},{2,4}})));
        Modelica.Blocks.Interfaces.RealInput BEox = inStream(bloodIn.BC[BCE.BEox])
          annotation (Placement(transformation(extent={{-8,-14},{2,-4}})));
      equation
       bloodIn.q = 0;
       bloodIn.BC = zeros(size(bloodIn.BC, 1)) "nothing goes out";

        connect(aBBB1.pH, totalO2Base1.pH) annotation (Line(points={{50,22.4},{48,22.4},
                {48,22},{54,22},{54,58},{-14,58},{-14,57},{-12,57}},
                                               color={0,0,127}));
        connect(aBBB1.Hb, Hb) annotation (Line(points={{18,22.4},{16,22.4},{16,22},{-100,
                22}},      color={0,0,127}));
        connect(aBBB1.pH, totalCO2Base1.pH) annotation (Line(points={{50,22.4},{48,22.4},
                {48,22},{54,22},{54,-30},{-2,-30},{-2,-31}},
                                            color={0,0,127}));
        connect(totalO2Base1.tO2, tO2) annotation (Line(points={{-30,61},{-100,
                61},{-100,62},{-100,60}},
                               color={0,0,127}));
        connect(totalO2Base1.pCO2, totalCO2Base1.pCO2) annotation (Line(points={{-20,49},
                {-18,49},{-18,20},{-18,-31}},          color={0,0,127}));
        connect(totalO2Base1.sO2, totalCO2Base1.sO2) annotation (Line(points={{-15,49},
                {-8,49},{-8,-36},{-8,-37},{-2,-37}}, color={0,0,127}));
        connect(totalCO2Base1.pO2, totalO2Base1.pO2) annotation (Line(points={{-2,-43},
                {-2,-43},{-2,-44},{62,-44},{62,67},{-11,67}}, color={0,0,127}));
        connect(totalCO2Base1.pCO2,aBBB1. pCO2) annotation (Line(points={{-18,-31},{-18,
                4},{-18,17.6},{18,17.6}},
                                    color={0,0,127}));
        connect(totalO2Base1.sO2, aBBB1.sO2)
          annotation (Line(points={{-15,49},{-8,49},{-8,14},{18,14},{18,12.8}},
                                                              color={0,0,127}));
        connect(aBBB1.P, P) annotation (Line(points={{18,6.4},{18,6},{18,6},{16,6},{16,
                7},{-3,7}},
              color={0,0,127}));
        connect(aBBB1.Alb, Alb)
          annotation (Line(points={{18,0},{18,0},{18,2},{18,2},{18,-1},{-3,-1}},
                                                                   color={0,0,127}));
        connect(BEox, aBBB1.BEox) annotation (Line(points={{-3,-9},{18,-9},{18,-6},{18,
                -6},{18,-6},{18,-6.4}},
                            color={0,0,127}));
        connect(totalCO2Base1.tCO2, tCO2) annotation (Line(points={{-20,-45},{-20,-45},
                {-20,-46},{-100,-46}}, color={0,0,127}));
        connect(Hb, totalO2Base1.Hb) annotation (Line(points={{-100,22},{-66,22},
                {-66,49},{-30,49}}, color={0,0,127}));
        connect(Hb, totalCO2Base1.Hb) annotation (Line(points={{-100,22},{-60,
                22},{-60,-41},{-20,-41}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics
              ={                                                     Text(
                extent={{-100,-120},{100,-20}},
                lineColor={0,0,127},
                fillColor={255,170,170},
                fillPattern=FillPattern.Solid,
                textString="%class"), Rectangle(extent={{-100,100},{100,-100}},
                  lineColor={0,255,0})}),                              Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end bloodABBMeasurementTestHack;
    end Thrash;
  end WBABB;
  annotation (uses(
      Physiomodel(version="0.2.29"),
      Physiolibrary(version="2.3.1"),
      Modelica(version="3.2.2")));
end FullBloodAcidBase;
