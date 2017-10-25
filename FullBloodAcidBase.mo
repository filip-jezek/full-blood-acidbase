within ;
package FullBloodAcidBase

  package FiggeFencl

  model FiggeFencl3
    "Base class for plasma acidbase after Figge and Fencl 3.0, missing inputs for SID, PCO2, total Pi and total albumin. Not meant to be run independently."

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
    input Real SID( displayUnit = "meq/l") "Strong ion difference. Normal value 39";
    input Real pCO2(displayUnit = "mmHg") "CO2 partial pressure. Normal value 40";
    input Real Pi( unit = "mmol/l") "Total phosphate. Normal value 1.15";
    input Real alb( unit= "g/dl") "Albumin concentration. Normal value 4.4";

    Real H( displayUnit="eq/l")= 10 ^ (-pH);
    Real pH(start = 10, unit="1");
    Real HO( displayUnit="eq/l") = kw / H;
    Real HCO3( displayUnit="eq/l") = Kc1 * pCO2 / H;
    Real CO3(  displayUnit="eq/l")= Kc2 * HCO3 / H;

    protected
    Real FNX = K1 * H^2 + 2 * K1 * K2 * H + 3 * K1 * K2 * K3;
    Real FNY = H^3 + K1 * H^2 + K1 * K2 * H + K1 * K2 * K3;
    Real FNZ = FNX / FNY;
    public
    Real P( displayUnit = "meq/l")= Pi * FNZ;
    Real Netcharge = SID + 1000 * (H - HO - HCO3 - CO3) - P;

    Real NB = 0.4 * (1 - (1 / (1 + 10 ^ (pH - 6.9))))
        "NB accounts for histidine pK shift due to the NB transition";

    // constant Real albuminResidues[:] = cat(1,{-1 /*cysteine */,-98/*glutamic acid*/,-18/*tyrosine*/,+24/*arginine */, /* lysine >>>*/ 2, 2, 2, 2, 1, 50} ,ones(16) /*histidine residues*/,/* amino terminus and carboxyl terminus*/{1, 1});
    protected
    Real albConversion = 1000 * 10 * alb / 66500;
    constant Real albuminResidues[:] = cat(1,{-1,              -98,                 -18,            +24,                              2, 2, 2, 2, 1, 50}, ones(16),                                                                {1, -1});
    // Real albuminPks[:] = {8.5 /* CYST*/,3.9 /* GLUT*/,11.7 /* TYR*/,12.5 /* ARG*/,/*LYS >>>*/5.8, 6.15, 7.51, 7.685, 7.86, 10.3,/*HIST>>>*/7.12 - NB, 7.22 - NB, 7.1 - NB, 7.49 - NB, 7.01 - NB, 7.31, 6.75, 6.36, 4.85, 5.76, 6.17, 6.73, 5.82, 5.1, 6.7, 6.2, 8/* amino terminus */,3.1 /*carboxyl terminus*/};
    Real albuminPks[:] = {8.5,          3.9,          11.7,         12.5,                    5.8, 6.15, 7.51, 7.685, 7.86, 10.3,           7.12 - NB, 7.22 - NB, 7.1 - NB, 7.49 - NB, 7.01 - NB, 7.31, 6.75, 6.36, 4.85, 5.76, 6.17, 6.73, 5.82, 5.1, 6.7, 6.2, 8,                    3.1};
    Real albChrg[n](each displayUnit = "meq/l") "charge of albumin per unit";
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
</html>",   revisions="<html>
<pre><font style=\"color: #006400; \">Filip Jezek, 2016</font></pre>
</html>"));
  end FiggeFencl3;

  model FiggeFencl3Detailed "Extension for investigation of detailed albumin balance"
    extends FiggeFencl.FiggeFencl3;
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
</html>",   revisions="<html>
<pre><font style=\"color: #006400; \">Filip Jezek, 2016</font></pre>
</html>"));
  end FiggeFencl3Detailed;

    model FiggeFenclNSID
      "Calculation of normal SID using plasma figge fencl"
      parameter Real pH0(displayUnit = "mmHg") = 7.4;
      input Real pCO20(displayUnit = "mmHg") = 40;
      input Real Pi0( unit = "mmol/l");//= 1.15;
      input Real alb0( unit="g/dl");//= 4.4;
      Real SID;// = figgeFencl3.SID;

      FiggeFencl.FiggeFencl3 figgeFencl3(
        pH=pH0,
        pCO2=pCO20,
        Pi=Pi0,
        alb=alb0,
        SID=SID) annotation (Placement(transformation(extent={{-58,0},{-38,20}})));
    end FiggeFenclNSID;

    model FiggeFencl3pCO2 "Change of pH while changing pCO2"
      extends Modelica.Icons.Example;

      FiggeFencl.FiggeFencl3Detailed figgeFencl3Base(
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

  model SimplePlasma
    "Base class for plasma acidbase after Figge and Fencl 3.0, missing inputs for SID, PCO2, total Pi and total albumin. Not meant to be run independently."

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
    input Real SID( displayUnit = "meq/l") "Strong ion difference. Normal value 39";
    input Real pCO2(displayUnit = "mmHg") "CO2 partial pressure. Normal value 40";
    input Real Pi( unit = "mmol/l") "Total phosphate. Normal value 1.15";
    input Real alb( unit= "g/dl") "Albumin concentration. Normal value 4.4";

    Real H( displayUnit="eq/l")= 10 ^ (-pH);
    Real pH(start = 10, unit="1");
    Real HO( displayUnit="eq/l") = kw / H;
    Real HCO3( displayUnit="eq/l") = Kc1 * pCO2 / H;
    Real CO3(  displayUnit="eq/l")= Kc2 * HCO3 / H;

    Real ZPi = (-1) - 10 ^ (pH - 6.87) / (1 + 10 ^ (pH - 6.87));
    Real ZAlb = (-10.7) - 16 * (10 ^ (pH - 7.42) / (1 + 10 ^ (pH - 7.42)));

    public
    Real P( displayUnit = "meq/l")= Pi * ZPi;
    Real Netcharge = SID + 1000 * (H - HO - HCO3 - CO3) - P;
    // constant Real albuminResidues[:] = cat(1,{-1 /*cysteine */,-98/*glutamic acid*/,-18/*tyrosine*/,+24/*arginine */, /* lysine >>>*/ 2, 2, 2, 2, 1, 50} ,ones(16) /*histidine residues*/,/* amino terminus and carboxyl terminus*/{1, 1});
    Real atch = ZAlb*1000 * 10 * alb / 66500
        "albumin total charge. Normal value -12.2678";

  equation
    Netcharge + atch = 0;

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
</html>",   revisions="<html>
<pre><font style=\"color: #006400; \">Filip Jezek, 2016</font></pre>
</html>"));
  end SimplePlasma;

    model FiggeFenclNSID_Simple
      "Calculation of normal SID using plasma figge fencl"
      parameter Real pH0(displayUnit = "mmHg") = 7.4;
      input Real pCO20(displayUnit = "mmHg") = 40;
      input Real Pi0( unit = "mmol/l");//= 1.15;
      input Real alb0( unit="g/dl");//= 4.4;
      Real SID;// = figgeFencl3.SID;

      FiggeFencl.SimplePlasma figgeFencl3(
        pH=pH0,
        pCO2=pCO20,
        Pi=Pi0,
        alb=alb0,
        SID=SID) annotation (Placement(transformation(extent={{-58,0},{-38,20}})));
    end FiggeFenclNSID_Simple;
  end FiggeFencl;

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
      Real d(displayUnit = "mEq/l") = BEox - 0.3*(1-sO2)
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

    model CombinedModel
      "Test combined model of Figge-fencl plasma and SA full hemoatocrite"

      replaceable FiggeFencl.FiggeFencl3Detailed plasma(
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
      FiggeFencl.FiggeFenclNSID normalPlasma(
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
    end CombinedModel;

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
          SAnomogram_formalization.SAoriginal plasmaSA(
            Hb=0,
            BEox=BE,
            pCO2=pCO2) "SA original with Hct 0 for comparison with FF"   annotation (Placement(transformation(extent={{60,20},{80,40}})));

          SAnomogram_formalization.SAoriginal normalBloodSA(
            Hb = Hb,
            BEox=BE,
            pCO2=pCO2) "SA original for comparison only"   annotation (Placement(transformation(extent={{60,20},{80,40}})));


        end ResultSetAtBE;

        model SetAtAlb

          CombinedModel full_blood_combo(
            BE=BE,
            pCO2=pCO2,
            alb=alb,
            Hb=Hb,
            Pi=Pi)
            annotation (Placement(transformation(extent={{20,20},{40,40}})));

          FiggeFencl.FiggeFencl3Detailed FF_plasma_only(
            SID=normalPlasma.SID + BE,
            pCO2=pCO2,
            Pi=Pi,
            alb=alb) "Plasma model acc to FiggeFencl for comparison only"
            annotation (Placement(transformation(extent={{60,-20},{80,0}})));

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
          FiggeFencl.FiggeFenclNSID normalPlasma(
            pH0=7.4,
            pCO20=40,
            alb0=alb,
            Pi0=Pi)
            annotation (Placement(transformation(extent={{-80,40},{-60,60}})));

        end SetAtAlb;

      end Auxiliary;

    end comparisson;

    model Wolf_full_blood "Implementation of the Wolf acidbase model, reduced to Plasma and Erythrocyte"
      Real a = mClP + HCO3P;
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
      Real MCl(unit="mmol") = m0ClE + m0ClP;
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
       Real Qpp[:] = {C_NaP,  + C_KP,  + 2 * C_CaP,  + 2 * C_MgP,  - C_ClP,  + ZPi * C_PiP,  + ZAlb * C_AlbP,  + ZimP};
      parameter Real X(unit="mEq")=0;
      //
      input Real BE = 0;
      Real BEe= (1 - 0.023 * 9) * (HCO3P - 24.4 + (2.3 * 9 + 7.7) * (pHP - 7.4));
    // Real SID = (1 - (1 - HCO3E / HCO3P) * fH * fB) * HCO3P + (1 - fH * fB) * (C_AlbP * (8 * pHP - 41) + C_PiP * (0.3 * pHP - 0.4)) + C_Hb * fB * (10.2 * pHP - 73.6) + C_DPG * fH * fB * (0.7 * pHP - 0.5);
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
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end Wolf_full_blood;

    model CombinedModelSimple
      "Test combined model of Figge-fencl plasma and SA full hemoatocrite"

      replaceable FiggeFencl.SimplePlasma plasma(
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
      FiggeFencl.FiggeFenclNSID_Simple normalPlasma(
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
    end CombinedModelSimple;
  end Full_Blood;

  package Figures
    extends Modelica.Icons.ExamplesPackage;

    model Figure1_BE_curve "BE curve for graphs with SA nomogram"
      extends Modelica.Icons.Example;
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

    end Figure1_BE_curve;

    model Figure_1
      extends Modelica.Icons.Example;
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
    end Figure_1;

    model Figure_2
      extends Modelica.Icons.Example;

      parameter Real Pi = 1.15;
      parameter Real alb = 4.4;
      parameter Real Hb = 15;
      Real BE = time*40 - 20;
      Real BEfixed = 0;
      Real pCO2 = time*60 + 20;
      Real pCO2fixed = 40;

      FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb fixedBE(  BE = BEfixed, pCO2 = pCO2, alb = alb, Hb = Hb, Pi = Pi)
        "Normal model set at fixed BE";
      FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb fixedpCO2(  BE = BE, pCO2 = pCO2fixed, alb = alb, Hb = Hb, Pi = Pi)
        "Normal model set at fixed pCO2";

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

    end Figure_2;

    model Figure_3
      "Dilution of FF and full blood models. The dilution of SA is made possible by incorporating the physicochemical model"
      extends Modelica.Icons.Example;
      parameter Real pCO2( displayUnit = "mmHg") = 40;
      Real BE( displayUnit = "meq/l") = ndp.SID*dilutionFactor - normalPlasma.SID "Shift in BE for diluted system.";
      Real newSID(displayUnit = "meq/l") = ndp.SID*dilutionFactor "SID to which the system is diluted to. Only for comparison.";
      // normal SID is 38.97
      parameter Real alb( unit = "g/dl")= 4.4;
      parameter Real Pi( unit = "mmol/l")= 1.15;
      Real dilutionFactor = 1.25 - time;
      Real Hb(  unit = "g/dl") = 15*dilutionFactor;
      Real HbPerCent = Hb/15*100;
      Full_Blood.comparisson.Auxiliary.SetAtAlb dilutedClosed(
        BE=BE,
        pCO2=pCO2*dilutionFactor,
        alb=alb*dilutionFactor,
        Hb=Hb,
        Pi=Pi*dilutionFactor)
        "Diluted closed system with diluted pCO2"
        annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
      Full_Blood.comparisson.Auxiliary.SetAtAlb dilutedOpenSystem(
        BE=BE,
        pCO2=pCO2,
        alb=alb*dilutionFactor,
        Hb=Hb,
        Pi=Pi*dilutionFactor) "Diluted open system, i.e. the pCO2 is held constant by the respiration"
        annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));
      Full_Blood.comparisson.Auxiliary.SetAtAlb nondilluted(
        BE=BE,
        pCO2=pCO2,
        alb=alb,
        Hb=Hb,
        Pi=Pi) "Non-dilluted system for comparison"
        annotation (Placement(transformation(extent={{0,-40},{20,-20}})));

      FiggeFencl.FiggeFenclNSID plasma(
        pH0=7.4,
        pCO20=pCO2*dilutionFactor,
        alb0=alb*dilutionFactor,
        Pi0=Pi*dilutionFactor);

    protected
      FiggeFencl.FiggeFenclNSID normalPlasma(
        pH0=7.4,
        pCO20=40,
        alb0=alb*dilutionFactor,
        Pi0=Pi*dilutionFactor)
        "normal plasma for establishing NSID in diluted blood";

      FiggeFencl.FiggeFenclNSID ndp(
        pH0=7.4,
        pCO20=40,
        alb0=alb,
        Pi0=Pi)
        "normal plasma for establishing NSID in reference, i.e. non-diluted blood";

    public
      FiggeFencl.FiggeFencl3 figgeFencl3Open(
        SID=ndp.SID*dilutionFactor,
        pCO2=pCO2,
        Pi=Pi*dilutionFactor,
        alb=alb*dilutionFactor)
        "Comparison to dilution of FiggeFencl3 model in open system, i.e. with constant pCO2"
        annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
      FiggeFencl.FiggeFencl3 figgeFencl3Closed(
        SID=ndp.SID*dilutionFactor,
        pCO2=pCO2*dilutionFactor,
        Pi=Pi*dilutionFactor,
        alb=alb*dilutionFactor)
        "Comparison to dilution of FiggeFencl3 model in closed system, i.e. with diluted pCO2"
        annotation (Placement(transformation(extent={{-80,0},{-60,20}})));

        // dilution graph in Dymola
    /* 
createPlot(id=6, position={46, 289, 586, 421}, x="dilutionFactor", y={"compensatedpCO2.FF_plasma_only.pH", "dilluted.FF_plasma_only.pH", 
"dilluted.full_blood_combo.pH", "compensatedpCO2.full_blood_combo.pH"}, range={0.5, 1.25, 7.1000000000000005, 7.500000000000001}, grid=true, legend=false, filename="dsres.mat", bottomTitleType=2, bottomTitle="Dilution factor", colors={{0,0,0}, {0,0,0}, {238,46,47}, {238,46,47}}, patterns={LinePattern.Solid, LinePattern.Dash, LinePattern.Dash, LinePattern.Solid}, thicknesses={0.5, 0.5, 0.5, 0.5});
*/

    end Figure_3;

    model Figure_5 "Effects of acute hypoalbuminemia"
      extends Modelica.Icons.Example;

      Real BEf = 0 "BE for fixed BE during Alb change";
      Real BEsidf = normalPlasma.SID - albPlasma.SID "BE for fixed SID during Alb change";

      parameter Real pCO2( displayUnit = "mmHg") = 40;
      parameter Real Hb( unit = "g/dl") = 15;
      parameter Real Pi( unit = "mmol/l") = 1.15;
      parameter Real alb0( unit = "g/dl") = 4.4;
      Real alb(  unit = "g/dl") = alb0/2 + alb0*time "Variable alb from (0.5 - 1.5)*alb0 during 1s";

      FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb fixedBE(  BE = BEf, pCO2 = pCO2, alb = alb, Hb = Hb, Pi = Pi)
        "Normal model set with fixed BE, that is Albumin is added (or removed) with according H+, resulting in HCO3 increase."
      annotation (Placement(transformation(extent={{20,20},{40,40}})));

      FullBloodAcidBase.Full_Blood.comparisson.Auxiliary.SetAtAlb fixedSID(  BE = BEsidf, pCO2 = pCO2, alb = alb, Hb = Hb, Pi = Pi)
        "Normal model set with fixed SID, that is Albumin is added (or removed) with dissociated cation (e.g. Na, K..)"
      annotation (Placement(transformation(extent={{20,20},{40,40}})));

    protected
      FiggeFencl.FiggeFenclNSID albPlasma(
        pH0=7.4,
        pCO20=40,
        alb0=alb,
        Pi0=Pi)
        "Normal SID for changed albumin level";

      FiggeFencl.FiggeFenclNSID normalPlasma(
        pH0=7.4,
        pCO20=40,
        alb0=alb0,
        Pi0=Pi)
        "Normal standard SID";

          // plot
          /*
createPlot(id=8, position={-39, 211, 865, 629}, x="alb", y={"fixedSID.FF_plasma_only.barGraphHCO3Alb1", "fixedSID.FF_plasma_only.barGraphHCO3Alb2",
 "fixedSID.FF_plasma_only.barGraphHCO3Alb3", "fixedBE.FF_plasma_only.barGraphHCO3Alb1",
 "fixedBE.FF_plasma_only.barGraphHCO3Alb2", "fixedBE.FF_plasma_only.barGraphHCO3Alb3",
 "fixedBE.FF_plasma_only.pH", "fixedSID.FF_plasma_only.pH"}, range={2.2, 6.6000000000000005, 0.0, 60.0}, grid=true, legend=false, filename="dsres.mat", leftTitleType=2, leftTitle="Charge [mEq/l]", bottomTitleType=2, bottomTitle="Albumin concentration [g/dl]", colors={{0,0,0}, {0,0,0}, {0,0,0}, {238,46,47}, {238,46,47}, {238,46,47}, {238,46,47}, 
{0,0,0}}, patterns={LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, LinePattern.Dot, 
LinePattern.Dot, LinePattern.Dot, LinePattern.Dash, LinePattern.Dash}, thicknesses={0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0}, range2={6.0, 7.6}, rightTitleType=2, rightTitle="pH", axes={1, 1, 1, 1, 1, 1, 2, 2});
*/
    end Figure_5;

    model Supplement_Figure1 "Comparison of SA nomogram formalisations"
        extends Modelica.Icons.Example;
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
      FullBloodAcidBase.Figures.Figure1_BE_curve BE_curve;

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
    end Supplement_Figure1;
  end Figures;

  package Tests
    extends Modelica.Icons.ExamplesPackage;

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

    model Test_Combo_Wolf
      Real pCO2 = 20 + time*60;
      Real BE = 0;
      parameter Real a = 0.65/0.94;
      Full_Blood.CombinedModel full_blood_combo(pCO2=pCO2, BE=BE)
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

  package Wolf
    package OriginalValues
      package Auxiliary
        record Erythrocyte
          import Modelica.SIunits.*;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));

          type cont = enumeration(
              Na,
              K,
              Cl,
              Hb,
              DPG,
              ATP,
              GSH,
              im,
              Pi) "Contents of erythrocyte";
          parameter Concentration Na=10;
          parameter Concentration K=99;
          parameter Concentration Cl0=53.8;
          Concentration Cl(start=Cl0);
          Real Cl_mass0(unit = "mol") = Cl0 / wf0 * water_volume0;
          Real Cl_mass(unit = "mol") = Cl / wf0 * water_volume;

          parameter Concentration Hb=5.3 "concentration of Hb tetramer";
          parameter Concentration DPG=4.4;
          parameter Concentration ATP=1.8;

          parameter Concentration GSH=2.2;
          parameter Concentration im=21.64 "adjusted value";
          // Concentration im=20.2 "Original value in article [1]";
          //   Concentration im;//=20.2;
          parameter Concentration Pi=0.67 "mono- and di- valent";

          Concentration volume_c[cont]={Na,K,Cl,Hb,DPG,ATP,GSH,im,Pi}
            "concentration in one liter";
          Concentration water_c[cont]=volume_c ./ wf0 .* water_volume0 ./ water_volume
            "Actual concentration recalculated to water fraction (initially 0.73) per one liter of volume";
        //   Real masses0[cont](each unit="mol") = volume_c ./ wf0 .* water_volume0
        //     "total mass of contents";

          Real ZHb=15.6 - 23*(10^(pH - 6.69)/(1 + 10^(pH - 6.69))) - 4*(10^(pH - 7.89)/(
              1 + 10^(pH - 7.89))) "For spO2=100, otherwise Zhb += 1.5*((1 - 0.75)/0.75)";
          Real ZDPG=(-3) - 1*(10^(pH - 7.56)/(1 + 10^(pH - 7.56))) - 1*(10^(pH - 7.32)/(
              1 + 10^(pH - 7.32)));
          Real ZATP=(-3) - 1*(10^(pH - 6.8)/(1 + 10^(pH - 6.8)));
          Real ZGSH=(-1) - 1*(10^(pH - 8.54)/(1 + 10^(pH - 8.54))) - 1*(10^(pH - 9.42)/(
              1 + 10^(pH - 9.42)));
          Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87))
            "this was not in the original implementatio by  machek!";
        //   Real Zim;//=-9.2 "Charge of ALL impermeable solutes";
          parameter Real Zim=-4.85 "Adjusted Charge of ALL impermeable solutes";

          //osmolality
          Real fiHb=1 + 0.115*water_c[cont.Hb] + 0.0256*water_c[cont.Hb]^2
            "Osmotic coefficient of Hemoglobin (eq from Raftos et al)";
          Real fi[cont]={0.93,0.93,0.93,fiHb,1,1,1,1,0.93} "Osmotic coefficients";
          Real Osm=sum(water_c .* fi) + 0.93*(HCO3 + CO3) "osmotic coefficient from Wolf model V3.49";
          //   //electric charge
          Real Z[cont]={1,1,-1,ZHb,ZDPG,ZATP,ZGSH,0,ZPi}
            "equivalent charges per mole";
          Real charge=(sum(water_c .* Z) - HCO3 - 2*CO3 + Zim)*water_volume;

          // parameter Fraction Hct0 = 0.44;
          // parameter Volume totalBlodVolume;
          parameter Volume water_volume0=0.44*0.73*1;
          //Hct * wf0 * Vblood;
          constant VolumeFraction wf0=0.73 "Fraction of water in the erythrocyte";

          Concentration HCO3=0.026*pCO2mmHg*10^(pH - 6.11) "HCO3 should be 13.6 acc to article";
          Concentration CO3=HCO3*10^(pH - 10.2);
          Real pH(start=7.22) = -log10(H);

          Volume water_volume(start=water_volume0);
          Real pCO2mmHg(unit="1");
          Real H(start=10^(-7.2));
        end Erythrocyte;

        record Plasma
          import Modelica.SIunits.*;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));

                type cont = enumeration(
              Na,
              K,
              Ca,
              Mg,
              Cl,
              Pi,
              Alb,
              im)
                "Contents of erythrocyte";

          parameter Concentration Na = 140;
          parameter Concentration K = 4.1;
          parameter Concentration Ca = 2.3;
          parameter Concentration Mg = 0.8;
          parameter Concentration Cl0 = 105;
          Concentration Cl( start = Cl0);
          Real Cl_mass0(unit = "mol") = Cl0 / wf0 * water_volume0;
          Real Cl_mass(unit = "mol") = Cl / wf0 * water_volume;
          parameter Concentration Pi = 1.2;
          parameter Concentration Alb = 0.65;
          // parameter Concentration im = 0; // original
          parameter Concentration im = 11.87; // adjusted

          Concentration volume_c[cont] = {Na, K, Ca, Mg, Cl, Pi, Alb, im} "concentration in one liter";
          Concentration water_c[cont]=volume_c / wf0 * water_volume0 / water_volume
            "Actual concentration recalculated to water fraction in one liter (initially 0.73)";
        //   Real masses0[cont](each unit="mol") = volume_c ./ wf0 .* water_volume0
        //     "total mass of contents";

          //charge on inpermeable solutes
          Real ZPi = (-1) - 10 ^ (pH - 6.87) / (1 + 10 ^ (pH - 6.87));
          Real ZAlb = (-10.7) - 16 * (10 ^ (pH - 7.42) / (1 + 10 ^ (pH - 7.42)));
        //   parameter Real Zim ;//= -5.3 "Charge of ALL impermeable solutes";//test
          // Real Zim = -5.3 "Charge of ALL impermeable solutes";//test
         parameter Real Zim = -10.99 "ADJUSTED Charge of ALL impermeable solutes";//test

          Real fi[cont]={0.93, 0.93, 1, 1, 0.93, 0.93, 1, 1}; // Alb has osmotic coefficient of 1???
          Real Osm=sum(water_c .* fi) + 0.93*(HCO3 + CO3);

          Real Z[cont]={1, 1, 2, 2, -1, ZPi, ZAlb, 0};
          Real Chrgs[cont] = water_c .* Z;
          Real charge=(sum(water_c .* Z) + Zim - HCO3 - 2*CO3);//*water_volume;

          parameter Volume water_volume0 = (1-0.44)*0.94*1;

          constant VolumeFraction wf0 = 0.94 "Fraction of water in the plasma";
          Concentration HCO3 = 0.0306 * pCO2mmHg * 10 ^ (pH - 6.11);
          Concentration CO3=HCO3*10^(pH - 10.2);
          Real pH(start=7.37) = -log10(H);

          Volume water_volume( start = water_volume0);
          Real pCO2mmHg(unit="1");
          Real H(start=10^(-7.37));
        end Plasma;
      end Auxiliary;

      model PE "Plasma - erythrocyte model"
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));

        FullBloodAcidBase.Wolf.OriginalValues.Auxiliary.Erythrocyte ery;
        FullBloodAcidBase.Wolf.OriginalValues.Auxiliary.Plasma pla;

      // total mass of Cl mobile ion
      Real totalCl = ery.Cl_mass0 + pla.Cl_mass0;
      Real totalWater = ery.water_volume0 + pla.water_volume0;
       Real pCO2mmHg = 40;
      equation
        // conservation of total mass of mobile ions
        ery.Cl_mass0 + pla.Cl_mass0 = ery.Cl_mass + pla.Cl_mass;

        // volume conservation
        ery.water_volume0 + pla.water_volume0 = ery.water_volume + pla.water_volume;

        // same osmolarities
        ery.Osm = pla.Osm;

        // zero charges
        ery.charge = 0;
        pla.charge = 0;

        // Donnan equilibrium
        ery.water_c[ery.cont.Cl] / pla.water_c[pla.cont.Cl] = pla.H / ery.H;

        // pco2
        ery.pCO2mmHg = pCO2mmHg;
        pla.pCO2mmHg = pCO2mmHg;

      end PE;

      package Tests
        model P
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(
                  preserveAspectRatio=false)));
          FullBloodAcidBase.Wolf.OriginalValues.Auxiliary.Plasma ery;

        // total mass of Cl mobile ion
        Real pCO2mmHg = 40;
        equation
          // volume conservation
          ery.water_volume0 = ery.water_volume;

          // same osmolarities
          ery.Osm = 285;
          ery.charge = 0;
          ery.Cl = ery.Cl0;

          // Donnan equilibrium
          ery.pH = 7.37;
          // pco2
          ery.pCO2mmHg = pCO2mmHg;

        end P;

        model E
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(
                  preserveAspectRatio=false)));
          FullBloodAcidBase.Wolf.OriginalValues.Auxiliary.Erythrocyte ery;

        // total mass of Cl mobile ion
        Real pCO2mmHg = 40;
        equation
          // volume conservation
          ery.water_volume0 = ery.water_volume;

          // same osmolarities
            ery.Osm = 285;
          ery.charge = 0;
          ery.Cl = ery.Cl0;

          // Donnan equilibrium
          ery.pH = 7.22;
          // pco2
          ery.pCO2mmHg = pCO2mmHg;

        end E;
      end Tests;
    end OriginalValues;

    package CurrentVersion
      package Auxiliary
        record Erythrocyte
          import Modelica.SIunits.*;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));

          type cont = enumeration(
              Na,
              K,
              Cl,
              Hb,
              DPG,
              ATP,
              GSH,
              im,
              Pi) "Contents of erythrocyte";
          parameter Concentration Na=13.69*wf0;
          parameter Concentration K=136*wf0;
          parameter Concentration Cl0=53.8;
          Concentration Cl(start=Cl0);
          Real Cl_mass0(unit="mol") = Cl0/wf0*Vew0;
          Real Cl_mass(unit="mol") = Cl/wf0*Vew;

          parameter Concentration Hb=5.3 "concentration of Hb tetramer";
           Concentration DPG=4.3/few*wf0 "the 4.3 value goes directly in";
          Concentration ATP=1.8/few*wf0;

           Concentration GSH=2.23/few*wf0;
          parameter Concentration im=16.03*wf0 "original value in model";
        //   parameter Concentration im=21.64 "adjusted value";
          // Concentration im=20.2 "Original value in article [1]";
          //   Concentration im;//=20.2;
          parameter Concentration Pi=0.67 "mono- and di- valent";
          parameter Real O2Sat(unit="1") = 0.75;
          Real Ve0ByVeFraction;
          Real fH;
          Real few;
          Real HbGperLiterBlood = Hb*Ve0ByVeFraction*fH*64.5;
          Concentration volume_c[cont]={Na,K,Cl,Hb,DPG,ATP,GSH,im,Pi}
            "concentration in one liter";
          Concentration water_c[cont]=volume_c ./ wf0 .* Vew0 ./ Vew
            "Actual concentration recalculated to water fraction (initially 0.73) per one liter of volume";
        //   Real masses0[cont](each unit="mol") = volume_c ./ wf0 .* water_volume0
        //     "total mass of contents";

          Real ZHb=15.6 - 23*(10^(pH - 6.69)/(1 + 10^(pH - 6.69))) - 4*(10^(pH - 7.89)/(
              1 + 10^(pH - 7.89))) + (1-O2Sat)*HbGperLiterBlood*0.012
                                                               "in article is Zhb += 1.5*((1 - 0.75)/0.75) in venous blood";
          Real ZDPG=(-3) - 1*(10^(pH - 7.56)/(1 + 10^(pH - 7.56))) - 1*(10^(pH - 7.32)/(
              1 + 10^(pH - 7.32)));
          Real ZATP=(-3) - 1*(10^(pH - 6.8)/(1 + 10^(pH - 6.8)));
          Real ZGSH=(-1) - 1*(10^(pH - 8.54)/(1 + 10^(pH - 8.54))) - 1*(10^(pH - 9.42)/(
              1 + 10^(pH - 9.42)));
          Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87))
            "this was not in the original implementatio by  machek!";
        //   Real Zim;//=-9.2 "Charge of ALL impermeable solutes";
        //   parameter Real Zim=-4.85 "Adjusted Charge of ALL impermeable solutes";
          parameter Real Zim=-8.0452 "Charge of ALL impermeable solutes acc to MODEL";

          //osmolality
          Real fiHb=1 + 0.115*water_c[cont.Hb] + 0.0256*water_c[cont.Hb]^2
            "Osmotic coefficient of Hemoglobin (eq from Raftos et al)";
          Real fi[cont]={0.93,0.93,0.93,fiHb,1,1,1,1,1} "Osmotic coefficients";
          Real Osm=sum(water_c .* fi) + 0.93*(HCO3 + CO3) + Lactate + permeableParticles "osmotic coefficient from Wolf model V3.49";
          Real OsmTest[:] =   water_c .* fi;
        //   //electric charge
          Real Z[cont]={1,1,-1,ZHb,ZDPG,ZATP,ZGSH,0,ZPi}
            "equivalent charges per mole";
          Real chargeTest[:] = water_c .* Z;
          Real charge=(sum(water_c .* Z) - HCO3 - 2*CO3 - Lactate + Zim*Vew0/Vew)*Vew;

          // parameter Fraction Hct0 = 0.44;
          // parameter Volume totalBlodVolume;
          Real Vew0;
          constant Real wf0(unit="1")=0.73 "Fraction of water in the erythrocyte";

          Concentration HCO3=0.026*pCO2mmHg*10^(pH - 6.103)/few "HCO3 should be 16,8 acc to the model";
          Concentration CO3=HCO3*10^(pH - 10.2) "CO3 is 0,016";
          Real pH(start=7.22) = -log10(H);

          Real Vew(unit="l", min=0, max=10, start = 1.43);
          Real pCO2mmHg(unit="1");
          Real H(start=10^(-7.2), min = 0, max = 1);
          Concentration Lactate "LACew =LACpw/rCl";
          parameter Concentration permeableParticles = 10.64 "glucose and urea concentration in PLasma water";

          Concentration HCO3ClDiff = HCO3/ water_c[cont.Cl]*few;
        end Erythrocyte;

        record Plasma
          import Modelica.SIunits.*;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));

          type cont = enumeration(
              Na,
              K,
              Ca,
              Mg,
              Cl,
              Pi,
              Alb,
              im,
              Lac) "Contents of erythrocyte";

          Concentration Na;
          Concentration K;
          Concentration Cl;

          Concentration Ca=2.4*wf0/fpw;
          parameter Concentration Mg=0.8;
          //parameter Concentration Cl0=105;
        //  Real Cl_mass0(unit="mol") = Cl0/wf0*water_volume0;
        //  Real Cl_mass(unit="mol") = Cl/wf0*water_volume;
          parameter Concentration Pi=1.2;
          Concentration Alb=AlbPwGPerL/66.5 "mmol/Lpw";
          Real AlbPwGPerL = 43 / Vp0ByVp "g/lpw";
          // parameter Concentration im = 0; // original
          parameter Concentration im=11.87;
          // adjusted
          parameter Concentration Lac=1.5;
          Concentration SO4pw=0.33/fpw "mmol/Lpw";
          Concentration volume_c[cont]={Na,K,Ca,Mg,Cl,Pi,Alb,im,Lac}
            "concentration in one liter";
          Concentration water_c[cont]=volume_c/fpw*Vp0ByVp
            "Actual concentration recalculated to water fraction in one liter (initially 0.96)";
          //   Real masses0[cont](each unit="mol") = volume_c ./ wf0 .* water_volume0
          //     "total mass of contents";
          Real CaPPBound=(3.98)/(3.98 + Hh)*water_c[cont.Ca] "1.257";
          Real MgPPBound=CaPPBound/2;
          Real ZCaBindPerAlb=CaPPBound/water_c[cont.Alb];
          Real ZClBindPerAlb=6.46/(6.46 + Hh)*6 + 4 "Anstey";
          Real CaIon=water_c[cont.Ca] - CaPPBound "mEq/Lpw";
          Real MgIon=water_c[cont.Mg] - MgPPBound "mEq/Lpw";
          //charge on inpermeable solutes
          Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87));
          Real ZFigge=(-10.65) - 16*(10^(pH - 7.418)/(1 + 10^(pH - 7.418)));
          Real ZAlbBnd = - ZClBindPerAlb + ZCaBindPerAlb + ZCaBindPerAlb/2;
          Real ZAlb=ZFigge + ZAlbBnd;
          Real Hh=H/fpw*1e8;

          //   parameter Real Zim ;//= -5.3 "Charge of ALL impermeable solutes";//test
          // Real Zim = -5.3 "Charge of ALL impermeable solutes";//test
          parameter Real Zim=-0.6349 "Charge of ALL impermeable solutes";

          Real fi[cont]={0.93,0.93,1,1,0.93,0.93,1,1,-9999};
          // Alb has osmotic coefficient of 1???
          Real CaOsm=water_c[cont.Ca] - CaPPBound*0.5;
          Real MgOsm=water_c[cont.Mg] - CaPPBound*0.5;
          Real OsmPart=SO4pw + CaOsm + MgOsm + sum(water_c[{cont.Na,cont.K,cont.Cl,cont.Lac}])
               + water_c[cont.Pi] + HCO3 + CO3;
          Real Osm=OsmPart*0.93 + permeableParticles;

          Real Z[cont]={1,1,2,2,-1,ZPi,ZAlb,0,1};
          Real Chrgs[cont]=water_c .* Z;

          Real SID=sum(water_c[{cont.Na,cont.K}]) - water_c[cont.Cl] + 2*sum(water_c[{cont.Ca,
              cont.Mg}]) - CaPPBound - MgPPBound - water_c[cont.Lac];
          Real charge=Zim + SID + water_c[cont.Alb]*ZAlb + water_c[cont.Pi]*ZPi - HCO3 -
              2*CO3 - SO4pw*2;
          //*water_volume;

          parameter Volume water_volume0=(1 - 0.44)*0.94*1;

          constant VolumeFraction wf0=0.94 "Fraction of water in the plasma";
          Concentration HCO3=0.0306*pCO2mmHg*10^(pH - 6.103)/fpw;
          Concentration CO3=HCO3*10^(pH - 10.2);
          Real pH(start=7.37) = -log10(H);

          Real Vp0ByVp;
          Real fpw "fraction plasma water - Vpw/Vp";
          Real pCO2mmHg(unit="1");
          Real H(start=10^(-7.37), min = 0, max = 1);
          Real Hw = H/fpw;
          parameter Concentration permeableParticles=10.64
            "glucose and urea concentration in PLasma water";
        end Plasma;

        record Isf
          "ISF is considered as water compartment, that is each substance concentration is meant as in water"
          import Modelica.SIunits.*;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));

          type cont = enumeration(
              Na,
              K,
              Ca,
              Mg,
              Cl,
              Pi,
              Alb,
              im,
              Lac) "Contents of erythrocyte";


        Concentration plasma_water_c[cont];

          Real rClpwis;
          // TODO own computation of ALB
          Concentration CaIon "mEq/l";
          Concentration MgIon                     "mEq/l";
          //Modelica.SIunits.Concentration CaPPBound "mEq/l";
          Concentration Cl;
          Real Vis;
          Real Vis0;
          parameter Real Alb(unit="g/l")= 43 "plasma alb conc g/l";
          Concentration SO4pw "mmol/Pw";

          Concentration albis = (Alb*0.33*Vis0)/(Vis - Vis0*0.25) "g/iw";
          Concentration albiw = albis/66.5 "mmol/liw Alb conc is fixed as 1/3 of 43, no matter the conc in plasma";

          Real CaPPBound= (3.98)/(3.98 + Hh)*plasma_water_c[cont.Ca] "= 1.290 Meq/Lis";
          Real CaPPBoundErr = 0 "Error in VisSim diagram, should be CaPPBound instead";
          Real ZCaBindPerAlb=CaPPBound/albiw;
          Real ZClBindPerAlb=6.46/(6.46 + Hh)*6 + 4 "Anstey";

          Real transf[cont] = {rClpwis, rClpwis, 0, 0, 0,(1/rClpwis)^(-ZPi), 0, 0, 1/rClpwis} "Only Na, K, Pi and Lac";

          Real SO4= (SO4pw*2)/(rClpwis^2)    "concentration in one liter mEq/lw";
          // Concentration water_c[cont]=transf .* plasma_water_c;
          Concentration water_c[cont]=transf .* plasma_water_c     "Actual concentration recalculated from plasma";

          Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87));
          Real ZFigge=(-10.65) - 16*(10^(pH - 7.418)/(1 + 10^(pH - 7.418)));
          Real ZAlbBnd = - ZClBindPerAlb + ZCaBindPerAlb + ZCaBindPerAlb/2;
          Real ZAlb=ZFigge+ ZAlbBnd;

          parameter Real Zim=3.057 "Charge of ALL impermeable solutes";

          Real OsmStep1=sum(water_c[{cont.Na,cont.K}]) + Cl;
          Real OsmStep2=OsmStep1 + (CaIon - CaPPBoundErr + MgIon - CaPPBoundErr/2)/2 + sum(water_c[
              {cont.Pi,cont.Lac}]) + HCO3 + CO3 + SO4/2;
          Real Osm = OsmStep2*0.93 + permeableParticles;

          Real SID=sum(water_c[{cont.Na,cont.K}]) - Cl + CaIon - CaPPBoundErr + (MgIon - 0.5
              *CaPPBoundErr) - water_c[cont.Lac];
          Real charge=Zim + SID + albiw*ZAlb + water_c[cont.Pi]*ZPi - HCO3 -
              2*CO3 - SO4;
          //*water_volume;

          Concentration HCO3=0.0326*pCO2mmHg*10^(pH - 6.103);
          Concentration CO3=HCO3*10^(pH - 10.2);
          Real pH(start=7.408) = -log10(H);

          Real pCO2mmHg(unit="1");
          Real H(start=10^(-7.408), min = 0, max = 1);
          Real Hh=H*1e8;

          parameter Concentration permeableParticles=10.64
            "glucose and urea concentration in PLasma water";


        // Transcapillary pressures
        // wattenpaugh J trauma Inuj 1998
          Real AlbPwGPerL;
          Real ca= AlbPwGPerL*0.1 "Albumin concentration";
          Real cp = AlbPwGPerL*0.1*1.637 "protein concentration";
          Real COPPpl =  PIEalb + PIEglob;
          Real PIEalb = ca / cp * ( cp*2.8+ 0.18*cp ^ 2 + 0.012*cp^3);
          Real PIEglob = (1 - ca / cp)*(0.9*cp + 0.12*cp^2 + 0.004*cp^3);
          Real cai = albis*0.1;
          Real COPalbis = cai*2.8 + 0.18*cai^2 + 0.012*cai^3;
        //   Real pd = (28.48 - 0.99*(COPPpl - COPalbis))/19.3;
          Real pd = (38.36 - 0.99*(COPPpl - COPalbis) + PVb)/19.3;
          Real OsmDiff = pd*2;
          Real Vbr = (Vb - Vb0)/Vb0 "VbRatio";
          Real PVb=if Vbr < -5 then -5*50 elseif Vbr > -5 and Vbr < 0 then Vbr*50
               elseif Vbr > 0 and Vbr < 100 then Vbr*46.5 else 100*46.5
            "In original -5..0 and 0..100";
          Real Vb;
          Real Vb0;

        end Isf;

        record Volumes
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));
          type PatientCat = enumeration(
              adultM,
              adultF,
              childM,
              childF);
          PatientCat pc=PatientCat.adultM;
          Real Fcell=if pc == PatientCat.adultM or pc == PatientCat.adultF then 0.88
               elseif pc == PatientCat.childM then 0.83 else 0.85;

                 // the rest is only adult male
           Real TBWm = 2.447 - 0.09156*age +0.3362* weight + 0.1074*height;
           //Real TBWf
           parameter Real LBM = 0.407*weight + 0.267*height - 19.2;
           // Real LBMm, LBMf;
           parameter Real age = 21;
           parameter Real weight = 70;
           parameter Real height = 175;

        parameter Real Ve0 = 0.048*LBM -0.72;
        parameter Real Vp0 = 0.051*LBM + 0.27;
        Real Vb0 = (Ve0 + Vp0);
        Real fH0 = Ve0 / (Ve0 + Vp0);
        Real fHa0 = fH0/Fcell;
        parameter Real Vis0 = 0.331*LBM - Vpw0;
        Real TBW = 2.447-0.09156*age+0.3362*weight+0.1074*height "Woodrow 2003 for male, there also exists female variant";
        Real Vc0 = TBW - Vis0 - Vpw0 - Vew0;


        parameter Real Vpw0= Vp0 - Vps;
        parameter Real Vps = 0.0603*Vp0;
        parameter Real Vew0 = Ve0*0.7317;
        parameter Real Ves = Ve0 - Vew0;
        Real Ve = Ves + Vew;
        Real Vb = Vpl + Ve;
        Real fH = Ve / Vb;

        Real Vw0 = Ve0 + Vp0 + Vis0 + Vc0;
        Real V = Vw0 + Vadd;
        parameter Real Vadd = 0;
        // Volumes
        Real Vpl = Vpw + Vps;
        Real fpw = Vpw / Vpl;
        Real few = Vew / Ve;
        Real Vpw =  Vadd + Vpw0 + Vis0 + Vew0 + Vc0 - Vis - Vew- Vc " = 2.93";
        Real Vp0ByVp = Vp0 / Vpl;
         Real Vew(start=Vew0) "= 1.43";
         Real Vis "= 15.59";
         Real Vc "= 22.9";
        // Real Vew = 1.43;
        // Real Vis = 15.59;
        // Real Vc = 22.9;

        end Volumes;

        record StrongIonMasses
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));

          Real Ve0;
          Real Vp0;
          Real Vis0;
          Real Vc0;
          // any ion addition MUST be added here as well;
          constant Real Napl0=140;
          constant Real Nais0=141.5;
          constant Real Nac0=12;
          Real MNa = Napl0*Vp0 + Nais0*Vis0 + Nac0*Vc0;
          Real MNa_IP;

          constant Real Kpl0=4.7;
          constant Real Kis0=4.75;
          constant Real Kc0=139;
          Real MK = Kpl0*Vp0 + Kis0*Vis0 + Kc0*Vc0;
          Real MK_IP;

          constant Real Cle0=53.5;
          constant Real Clpl0=104;
          constant Real Clis0=116.5;
          constant Real Clc0=4;
          Real MCl = Cle0*Ve0 + Clpl0*Vp0 + Clis0*Vis0 + Clc0*Vc0;


        end StrongIonMasses;

        record Cell
            import Modelica.SIunits.*;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));

          Real Kis;
          Real Nais;
          Real Clis;
          Real Vc0;
          Real Vc;
          Real His;
          Real pCO2mmHg(unit="1");


          Real Clc;

          Real rClisce = Clis/Clc;
          Real Kc = Kis*rClisce;
          Real Nac = Nais*rClisce*0.0029;

          Concentration HCO3=0.029*pCO2mmHg*10^(pH - 6.103);
          Real pH(start=7.408) = -log10(H);

          Real H(start=10^(-6.938), min = 0, max = 1) = His*0.1*rClisce "Why His*0.1? Is it the water content?";

          Real charge = -HCO3 + Kc - Clc +Nac + Zim*Im;
          Real Zim = -1.1728*10^(pH - 5.5)/(1 + 10^(pH - 5.5));
          Real Im = Vc0/Vc * 123.045;

          Real OsmPart = HCO3 + Kc+ Nac + Clc;
          Real Osm = Im + permeableParticles + 0.93*OsmPart;

          parameter Concentration permeableParticles=10.64
            "glucose and urea concentration in PLasma water";

        end Cell;
      end Auxiliary;

      model PE "Plasma - erythrocyte model"
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));

        FullBloodAcidBase.Wolf.CurrentVersion.Auxiliary.Erythrocyte ery;
        FullBloodAcidBase.Wolf.CurrentVersion.Auxiliary.Plasma pla;
        FullBloodAcidBase.Wolf.CurrentVersion.Auxiliary.Volumes vols;
        FullBloodAcidBase.Wolf.CurrentVersion.Auxiliary.StrongIonMasses sim;
        // total mass of Cl mobile ion
        Real pCO2mmHg=20 + time*40;
        Real rClpw_is=pla.water_c[pla.cont.Cl]/Clis "0.949290";
        constant Real MNac=276.96;
        constant Real MKc=3179.97;
        constant Real MClc=90.709;
        //constant Real MCle = 102.640;
        Real MCle=ery.water_c[ery.cont.Cl]*vols.Vew;
        constant Real Clis=116.79118;
        Real MCl_IP=sim.MCl - MClc - MCle;
        //Real test = sim.MNa_IP /(vols.Vis * rClpw_is + vols.Vpw);
        parameter Real Hct=0.44;
        parameter Real O2s=0.75;
        Real rClep=ery.water_c[ery.cont.Cl]/pla.water_c[pla.cont.Cl];
        Real rHep=pla.Hw/ery.H;
        Real Clpla( start = 110.8687, min = 0, max = 1000) = pla.water_c[pla.cont.Cl];
        Real Clery( start = 71.7743, min = 0, max = 1000) = ery.water_c[ery.cont.Cl];
      equation
        // vols.Vew = 1.43;
        vols.Vis = 15.5919;
        vols.Vc = 22.9;
        ery.Lactate = pla.water_c[pla.cont.Lac]*rClep;//1.035;
        //LACew =LACpw/rCl;
        //  ery.water_c[ery.cont.Cl] / 110.9 = 4.119084e-8 / ery.H "Clpla , Hpl at standard V3.49";
      //  ery.water_c[ery.cont.Cl]/pla.water_c[pla.cont.Cl] = pla.H/ery.H    "Clpla , Hpl at standard V3.49";

        rClep = rHep;
      pla.Osm = ery.Osm;
      //  pla.water_c[pla.cont.Cl] = 110.8687;
      //  ery.water_c[ery.cont.Cl] = 71.7743;



         ery.charge = 0;
         pla.charge = 0;
      //   pla.pH = 7.41225;
      //   ery.pH = 7.196;

        sim.Ve0 = vols.Ve0;
        sim.Vp0 = vols.Vp0;
        sim.Vis0 = vols.Vis0;
        sim.Vc0 = vols.Vc0;
        sim.MNa_IP = sim.MNa - MNac;
        sim.MK_IP = sim.MK - MKc;

        vols.Ve0/vols.Ve = ery.Ve0ByVeFraction;
        ery.fH = vols.fH;
        ery.few = vols.few;

        ery.Vew0 = vols.Vew0;
        ery.Vew = vols.Vew;
        pla.Vp0ByVp = vols.Vp0ByVp;
        pla.fpw = vols.fpw;


        // pco2
        pla.pCO2mmHg = pCO2mmHg;
        ery.pCO2mmHg = pCO2mmHg;


        pla.water_c[pla.cont.Na] = sim.MNa_IP/(vols.Vis*rClpw_is + vols.Vpw);
        pla.water_c[pla.cont.K] = sim.MK_IP/(vols.Vis*rClpw_is + vols.Vpw);
        Clis = (MCl_IP - vols.Vpw*pla.water_c[pla.cont.Cl])/vols.Vis;
        // pla.volume_c[pla.cont.Cl] = 110.86;

        // same osmolarities
        //   pla.Osm = 286.99;
        //   pla.charge = 0;



        // Donnan equilibrium




      end PE;

      package Tests
        model E
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(
                  preserveAspectRatio=false)));
          FullBloodAcidBase.Wolf.CurrentVersion.Auxiliary.Erythrocyte ery;
          FullBloodAcidBase.Wolf.CurrentVersion.Auxiliary.Volumes vols;
        // total mass of Cl mobile ion
        Real pCO2mmHg = 40;
        // output: osmolarity and Cl
        parameter Real Hct = 0.44;
        parameter Real O2s = 0.75;

        equation
          vols.Vew = 1.43;
          vols.Vis = 15.59;
          vols.Vc = 22.9;

          vols.Ve0/vols.Ve = ery.Ve0ByVeFraction;
          ery.fH = vols.fH;
          ery.few = vols.few;
          ery.Lactate = 1.035;//LACew =LACpw/rCl;
          ery.Vew0 = vols.Vew0;
          ery.Vew = vols.Vew;
        //  ery.water_c[ery.cont.Cl] / 110.9 = 4.119084e-8 / ery.H "Clpla , Hpl at standard V3.49";
          ery.water_c[ery.cont.Cl] = 71.77;
         ery.H = 6.362e-8;
          // ery.
        //   // same osmolarities
           //ery.Osm = 286.99;
         //  ery.charge = 0;
        //   ery.Cl = ery.Cl0;
           //ery.Cl = 71.77 * ery.wf0;

          // Donnan equilibrium
          // pco2
          ery.pCO2mmHg = pCO2mmHg;

        end E;

        model P
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(
                  preserveAspectRatio=false)));
          FullBloodAcidBase.Wolf.CurrentVersion.Auxiliary.Plasma pla;
          FullBloodAcidBase.Wolf.CurrentVersion.Auxiliary.Volumes vols;
          FullBloodAcidBase.Wolf.CurrentVersion.Auxiliary.StrongIonMasses sim;
        // total mass of Cl mobile ion
        Real pCO2mmHg = 40;
        // Real rClpw_is = pla.volume_c[pla.cont.Cl]/Clis;
        Real rClpw_is = 0.949290;
        constant Real MNac = 276.96;
        constant Real MKc = 3179.97;
        constant Real MClc = 90.709;
        constant Real MCle = 102.640;
        constant Real Clis = 116.79118;
        Real MCl_IP = sim.MCl - MClc - MCle;
        //Real test = sim.MNa_IP /(vols.Vis * rClpw_is + vols.Vpw);
        equation
          vols.Vew = 1.43;
          vols.Vis = 15.5919;
          vols.Vc = 22.9;

          sim.Ve0 = vols.Ve0;
          sim.Vp0 = vols.Vp0;
          sim.Vis0 = vols.Vis0;
          sim.Vc0 = vols.Vc0;
          sim.MNa_IP = sim.MNa - MNac;
          sim.MK_IP = sim.MK - MKc;


          pla.water_c[pla.cont.Na] = sim.MNa_IP / (vols.Vis * rClpw_is + vols.Vpw);
          pla.water_c[pla.cont.K] = sim.MK_IP / (vols.Vis * rClpw_is + vols.Vpw);
          Clis = (MCl_IP - vols.Vpw*pla.water_c[pla.cont.Cl])/ vols.Vis;
         // pla.volume_c[pla.cont.Cl] = 110.86;

          // same osmolarities
        //   pla.Osm = 286.99;
        //   pla.charge = 0;


          pla.Vp0ByVp = vols.Vp0ByVp;

          // Donnan equilibrium
          pla.pH = 7.41225;
          pla.fpw = vols.fpw;
          // pco2
          pla.pCO2mmHg = pCO2mmHg;

        end P;

        model I
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));
          FullBloodAcidBase.Wolf.CurrentVersion.Auxiliary.Isf  isf;

          Real rClpwis = 0.94929;
        equation

          isf.plasma_water_c={149.031, 4.962, 2.55814, 0, 0, 1.279073, 0, 0, 1.598841};

          isf.rClpwis = rClpwis;
          // TODO own computation of ALB
          isf.CaIon = 3.859*rClpwis^2;
          isf.MgIon = 1.0768*rClpwis^2;
          isf.Cl = 116.791;
          isf.Vis = 15.591;
          isf.Vis0 = 15.602;
          isf.SO4pw = 0.702420/2; // The input here is in mmol/lpw to match plaswma variable, not in mEq/lpw as in the vissim model
          isf.pCO2mmHg = 40;
          isf.pH = 7.4078;
          isf.AlbPwGPerL = 43.0654;
          isf.Vb = 5.08026;
          isf.Vb0 = 5.095485;
        end I;

        model C
          import Modelica.SIunits.*;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                coordinateSystem(preserveAspectRatio=false)));
                FullBloodAcidBase.Wolf.CurrentVersion.Auxiliary.Cell c;

        equation
          c.Kis = 4.7106;
          c.Nais = 141.474;
          c.Clis = 116.7911;
          c.Vc0 = 22.87176;
          c.Vc = 22.897789;
          c.His = 3.91020e-8;
          c.pCO2mmHg = 40;

          c.Clc = 3.961488;


        end C;
      end Tests;
      annotation (Documentation(info="<html>
<p>References:</p>
<p><br>Raftos, J. E., Bulliman, B. T., &amp; Kuchel, P. W. (1990). Evaluation of an electrochemical model of erythrocyte pH buffering using 31P nuclear magnetic resonance data. The Journal of General Physiology, 95(6), 1183&ndash;1204.</p>
<p>Wolf, M. B. (2013). Whole body acid-base and fluid-electrolyte balance: a mathematical model. American Journal of Physiology. Renal Physiology, 305(8), F1118&ndash;31.</p>
<p>Wolf, M. B. (2015). Comprehensive diagnosis of whole-body acid-base and fluid-electrolyte disorders using a mathematical model and whole-body base excess. Journal of Clinical Monitoring and Computing, 29(4), 475&ndash;490.</p>
<p>Wolf, M. B., &amp; DeLand, E. C. (2011). A comprehensive, computer-model-based approach for diagnosis and treatment of complex acid&ndash;base disorders in critically-ill patients. Journal of Clinical Monitoring and Computing, 25(6), 353&ndash;364.</p>
<p>Wolf, M. B., &amp; Deland, E. C. (2011). A mathematical model of blood-interstitial acid-base balance: application to dilution acidosis and acid-base status. Journal of Applied Physiology, 110(4), 988&ndash;1002.</p>
</html>"));
    end CurrentVersion;
  end Wolf;
  annotation (uses(
      Physiomodel(version="0.2.29"),
      Physiolibrary(version="2.3.1"),
      Modelica(version="3.2.2")));
end FullBloodAcidBase;
