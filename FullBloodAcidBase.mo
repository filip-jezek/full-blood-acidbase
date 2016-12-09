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
  Real albHAPlus = albConversion*sum(albPositivePart) "A0 + H+ = HA+";
  Real albAMinus = -albConversion*sum(albNegativePart) "A- + H+ = HA0";
  Real albA0 = (ATotPlus - albHAPlus) "A0 + H+ = HA+";
  Real albHA0 = (ATotMinus - albAMinus) "HA0 = A- + H+";
  Real ATotPlus = albConversion*sum(albTotalPlusPart)
      "Part of albumin, which could be positive";
  Real ATotMinus = -albConversion*sum(albTotalMinusPart)
      "Part of albumin, which could be negative.";

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
     input Real BEox(unit = "mEq/l") = 0;
     input Real sO2(unit = "1") = 1;
     input Real pCO2(unit = "mmHg") = 40;

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
       Real ctHb = Hct*33.34;
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

      FiggeFencl3Detailed                         plasma(
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
      parameter Real k = 1
        "pisvajcova konstanta - kolik mHCO3 zvostane v erytrocytu? Nebo jina nepresnost? Ideal 0.77";

      Real BEp( unit = "meq/l") = BE - mHCO3/(1-Hct);
      Real BEe( unit = "meq/l")= BE + k*mHCO3/Hct;

      Real mHCO3;
      //Real pHFullBlood;

      Real SID;
      input Real BE = 0;
    //  Real BE = time*40 - 20;
      input Real pCO2 = 40;
    //  Real logpco2 = log10(pCO2);
    //  input Real pCO2 = time*40 + 20;
      input Real Pi = 1.15;
      input Real alb = 4.4;
      output Real pH = fullErythrocyte.pH;
    //  Real bdHb = der(mHCO3)/der(plasma.pH);
    protected
      FiggeFenclNSID normalPlasma(
        pH0=7.4,
        pCO20=40,
        alb0=alb,
        Pi0=Pi)
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

            FullBloodAcidBase.SAnomogram_formalization.SAVanSlyke sAVanSlyke(  pCO2 = pCO2, BEox = BE, Hct = Hb/33.34, sO2 = 1, Alb= alb)
            annotation (Placement(transformation(extent={{-14,40},{6,60}})));
            input Real BE, pCO2;
            input Real df = 1 "dilution factor";
            input Real alb, Pi;
            input Real Hb;
        protected
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
  end Full_Blood;

  package Figures

    model Figure1
      Real pCO2 = time*40 + 20;
      SAnomogram_formalization.SA_comparison_pCO2 sA_comparison_pCO2_BE_10(BEox=
           -15, pCO2=pCO2)
        annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
      SAnomogram_formalization.SA_comparison_pCO2 sA_comparison_pCO2_BE0(BEox=0,
          pCO2=pCO2)
        annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
      SAnomogram_formalization.SA_comparison_pCO2 sA_comparison_pCO2_BE10(BEox=
            15, pCO2=pCO2)
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
    end Figure1;

    model Figure1_3_4_BE_curve
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

    end Figure1_3_4_BE_curve;

    model figure2A
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
      parameter Real safe_hct = 1;
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
    end figure2A;

    model figure2B
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
    end figure2B;

    model Figures3_4
      Real logpCO2 = log10(pCO2);
      Real pCO2 = time*40 + 20;
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

        /* FIGURE 3: Full blood combination compared to full blood SA, Figge of plasma and SA of Hct 0 (thus plasma as well)
    createPlot(id=3, position={70, 340, 483, 300}, x="pCO2", y={"resultSetAtBE_10.plasmaSA.pH", "resultSetAtBE_10.Normal.FF_plasma_only.pH", 
"resultSetAtBE_10.Normal.full_blood_combo.pH", "resultSetAtBE_10.normalBloodSA.pH",
 "resultSetAtBE0.plasmaSA.pH", "resultSetAtBE0.Normal.FF_plasma_only.pH", 
"resultSetAtBE0.Normal.full_blood_combo.pH", "resultSetAtBE0.normalBloodSA.pH",
 "resultSetAtBE10.plasmaSA.pH", "resultSetAtBE10.Normal.FF_plasma_only.pH", 
"resultSetAtBE10.Normal.full_blood_combo.pH", "resultSetAtBE10.normalBloodSA.pH"}, range={20.0, 60.0, 7.1000000000000005, 7.800000000000001}, autoscale=false, grid=true, legend=false, filename="BErange_alb_Hb5.mat", logX=true, leftTitleType=0, bottomTitleType=0, colors={{0,0,0}, {0,0,0}, {238,46,47}, {0,140,72}, {0,0,0}, {0,0,0}, {238,46,47}, 
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
 "resultSetAtBE_10.LowAlb.full_blood_combo.pH"}, range={20.0, 60.0, 7.1000000000000005, 7.800000000000001}, autoscale=false, grid=true, legend=false, filename="BErange_alb_Hb5.mat", logX=true, leftTitleType=0, bottomTitleType=0, colors={{238,46,47}, {28,108,200}, {0,0,0}, {28,108,200}, {0,0,0}, {238,46,47}, 
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
 "resultSetAtBE_10.LowAlb.full_blood_combo.pH", "resultSetAtBE_10.LowAlb.sAVanSlyke.pH"}, range={20.0, 60.0, 7.1000000000000005, 7.800000000000001}, autoscale=false, grid=true, legend=false, filename="dsres.mat", logX=true, leftTitleType=0, bottomTitleType=0, colors={{28,108,200}, {0,0,0}, {238,46,47}, {28,108,200}, {0,0,0}, {238,46,47}, 
{238,46,47}, {0,0,0}, {28,108,200}, {238,46,47}, {0,0,0}, {28,108,200}, 
{0,0,0}, {238,46,47}, {0,0,0}, {28,108,200}, {238,46,47}, {28,108,200}}, patterns={LinePattern.Dash, LinePattern.Dash, LinePattern.Dash, LinePattern.Solid, 
LinePattern.Solid, LinePattern.Solid, LinePattern.Dash, LinePattern.Dash, 
LinePattern.Dash, LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Dash, LinePattern.Solid, LinePattern.Solid, LinePattern.Solid, 
LinePattern.Dash, LinePattern.Dash}, thicknesses={0.25, 0.25, 0.5, 0.25, 0.25, 0.5, 0.5, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.5, 
0.25, 0.25, 0.5, 0.25}, rightTitleType=0);

*/
    end Figures3_4;

    model Figure5
      parameter Real Pi = 1.15;
      parameter Real alb = 4.4;
      parameter Real Hb = 15;
      Real BE = time*40 - 20;
      Real BEfixed = 0;
      Real pCO2 = time*40 + 20;
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
      Real BE= plasma.SID - normalPlasma.SID;
      Real newSID = ndp.SID*dilutionFactor;
      // mnoram normal SID is 38.97
      parameter Real alb = 4.4;
      parameter Real Pi = 1.15;
      Real dilutionFactor = 0.5 + time;
      Real Hb = 15*dilutionFactor;
      Full_Blood.comparisson.Auxiliary.SetAtAlb dilluted(
        df=dilutionFactor,
        BE=BE,
        pCO2=pCO2*dilutionFactor,
        alb=alb*dilutionFactor,
        Hb=Hb,
        Pi=Pi*dilutionFactor) "Dilluted"
        annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
      Full_Blood.comparisson.Auxiliary.SetAtAlb compensatedpCO2(
        df=dilutionFactor,
        BE=BE,
        pCO2=pCO2,
        alb=alb*dilutionFactor,
        Hb=Hb,
        Pi=Pi*dilutionFactor) "Dilluted"
        annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
      Full_Blood.comparisson.Auxiliary.SetAtAlb nondilluted(
        df=dilutionFactor,
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
  end Tests;

  package AlbuminBorderFlux

    model AlbuminBalance

    Physiolibrary.Chemical.Components.Diffusion UT_Capillary(Conductance(
          displayUnit="l/day") = 3.9351851851852e-09)
      annotation (Placement(transformation(extent={{-2,70},{10,82}})));
    Physiolibrary.Chemical.Components.Diffusion MT_Capillary(Conductance(
          displayUnit="l/day") = 7.4074074074074e-09)
      annotation (Placement(transformation(extent={{-2,38},{10,50}})));
    Physiolibrary.Chemical.Components.Diffusion LT_Capillary(Conductance(
          displayUnit="l/day") = 1.1805555555556e-08)
      annotation (Placement(transformation(extent={{0,12},{12,24}})));
    Physiolibrary.Chemical.Sources.UnlimitedSolutePump Transfusion(
          useSoluteFlowInput=false, SoluteFlow=0)
      annotation (Placement(transformation(extent={{20,-38},{0,-18}})));
    Physiolibrary.Chemical.Components.Stream UT_Lymph(useSolutionFlowInput=
            false, SolutionFlow=5.5333333333333e-09)
      annotation (Placement(transformation(extent={{10,66},{0,56}})));
    Physiolibrary.Chemical.Components.Stream MT_Lymph(useSolutionFlowInput=
            false, SolutionFlow=1.315e-08)
      annotation (Placement(transformation(extent={{10,34},{0,24}})));
    Physiolibrary.Chemical.Components.Stream LT_Lymph(useSolutionFlowInput=
            false, SolutionFlow=1.5933333333333e-08)
      annotation (Placement(transformation(extent={{10,8},{0,-2}})));
      Physiolibrary.Chemical.Components.Substance plasma(
        stateName="PlasmaProtein.Mass",
        useNormalizedVolume=false,
      solute_start=0.00437)
        annotation (Placement(transformation(extent={{-76,24},{-56,44}})));
      Physiolibrary.Chemical.Components.Substance UpperTorso(
        stateName="UT_InterstitialProtein.Mass",
        useNormalizedVolume=false,
        solute_start=0.00122)
        annotation (Placement(transformation(extent={{78,56},{58,76}})));
      Physiolibrary.Chemical.Components.Substance MiddleTorso(
        stateName="MT_InterstitialProtein.Mass",
        useNormalizedVolume=false,
      solute_start=0.00299)
        annotation (Placement(transformation(extent={{78,26},{58,46}})));
      Physiolibrary.Chemical.Components.Substance LowerTorso(
        stateName="LT_InterstitialProtein.Mass",
        useNormalizedVolume=false,
      solute_start=0.0018)
        annotation (Placement(transformation(extent={{78,-2},{58,18}})));
      AlbuminSynthesis                               synthesis(
          UseSythesisFactorInput=false, SynthesisBasic=1.6666666666667e-07)
        annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));
      Degradation                                      degradation(
        DegradationBasic=1.6666666666667e-07)
        annotation (Placement(transformation(extent={{0,-60},{20,-40}})));
    Physiolibrary.Chemical.Components.Diffusion GlomerulusProtein_Perm(
        Conductance=(0)*(1e-6)/60)
      annotation (Placement(transformation(extent={{0,-24},{20,-4}})));
      Physiolibrary.Chemical.Components.Substance Bladder(
        stateName="BladderProtein.Mass",
        useNormalizedVolume=false,
      solute_start=1e-15)
        annotation (Placement(transformation(extent={{78,-24},{58,-4}})));
    Physiolibrary.Chemical.Sensors.ConcentrationMeasure concentrationMeasure1
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-94,16})));
      ProteinDivision proteinDivision
        annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
      Physiolibrary.Chemical.Components.Clearance clearance(
          useSolutionFlowInput=true)
        annotation (Placement(transformation(extent={{-22,78},{-2,98}})));
      Modelica.Blocks.Sources.Pulse pulse(
        width=100,
        period(displayUnit="h") = 3600,
        nperiod=1,
        amplitude=1e-6,
        startTime(displayUnit="h") = 36000)
        annotation (Placement(transformation(extent={{22,82},{6,98}})));
      Physiolibrary.Types.Constants.VolumeConst volume(k=0.006063)
        annotation (Placement(transformation(extent={{98,36},{90,44}})));
      Physiolibrary.Types.Constants.VolumeConst volume1(k=0.00185)
        annotation (Placement(transformation(extent={{98,66},{90,74}})));
      Physiolibrary.Types.Constants.VolumeConst volume2(k=0.00247)
        annotation (Placement(transformation(extent={{98,8},{90,16}})));
      Physiolibrary.Types.Constants.VolumeConst volume3(k=0.0003)
        annotation (Placement(transformation(extent={{98,-14},{90,-6}})));
      Physiolibrary.Types.Constants.VolumeConst volume4(k=0.002807) annotation
        (Placement(transformation(
            extent={{-4,-4},{4,4}},
            rotation=0,
            origin={-74,52})));
      Physiolibrary.Types.Constants.pHConst pH(k=7.4)
        annotation (Placement(transformation(extent={{20,-84},{28,-76}})));
      Physiolibrary.Chemical.Sensors.MolarFlowMeasure molarFlowMeasure
        annotation (Placement(transformation(extent={{-54,-60},{-34,-40}})));
      Physiolibrary.Chemical.Sensors.MolarFlowMeasure molarFlowMeasure1
        annotation (Placement(transformation(extent={{-22,-60},{-2,-40}})));
      Modelica.Blocks.Math.Add add(k1=-1)
        annotation (Placement(transformation(extent={{0,-70},{10,-60}})));
      ProteinDivision proteinDivision1
        annotation (Placement(transformation(extent={{20,-70},{30,-60}})));
      ProteinCharge proteinCharge
        annotation (Placement(transformation(extent={{40,-80},{60,-60}})));
      AcidBaseBuffers acidBaseBuffers
        annotation (Placement(transformation(extent={{80,-80},{100,-60}})));
    equation
      connect(UT_Capillary.q_out,UpperTorso. q_out) annotation (Line(
          points={{10,76},{18,76},{18,66},{68,66}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,UT_Capillary. q_in) annotation (Line(
          points={{-66,34},{-26,34},{-26,76},{-2,76}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,UT_Lymph. q_out) annotation (Line(
          points={{-66,34},{-26,34},{-26,60},{0,60},{0,61}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,MT_Capillary. q_in) annotation (Line(
          points={{-66,34},{-26,34},{-26,42},{-14,42},{-14,44},{-2,44}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(MT_Capillary.q_out,MiddleTorso. q_out) annotation (Line(
          points={{10,44},{16,44},{16,36},{68,36}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(MT_Lymph.q_in,MiddleTorso. q_out) annotation (Line(
          points={{10,29},{16,29},{16,36},{68,36}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,MT_Lymph. q_out) annotation (Line(
          points={{-66,34},{-26,34},{-26,28},{0,28},{0,29}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,LT_Capillary. q_in) annotation (Line(
          points={{-66,34},{-26,34},{-26,18},{0,18}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,LT_Lymph. q_out) annotation (Line(
          points={{-66,34},{-26,34},{-26,2},{0,2},{0,3}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,GlomerulusProtein_Perm. q_in) annotation (Line(
          points={{-66,34},{-26,34},{-26,-14},{0,-14}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,Transfusion. q_out) annotation (Line(
          points={{-66,34},{-26,34},{-26,-28},{0,-28}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(LT_Capillary.q_out,LowerTorso. q_out) annotation (Line(
          points={{12,18},{18,18},{18,8},{68,8}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(LT_Lymph.q_in,LowerTorso. q_out) annotation (Line(
          points={{10,3},{18,3},{18,8},{68,8}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(GlomerulusProtein_Perm.q_out,Bladder. q_out) annotation (Line(
          points={{20,-14},{68,-14}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
      connect(plasma.q_out,concentrationMeasure1. q_in) annotation (Line(
          points={{-66,34},{-94,34},{-94,16}},
          color={107,45,134},
          thickness=1,
          smooth=Smooth.None));
    connect(concentrationMeasure1.concentration,proteinDivision. totalProteins)
      annotation (Line(
        points={{-94,24},{-94,10},{-80,10}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(UT_Lymph.q_in,UpperTorso. q_out) annotation (Line(
        points={{10,61},{18,61},{18,66},{68,66}},
        color={107,45,134},
        thickness=1,
        smooth=Smooth.None));
      connect(pulse.y, clearance.solutionFlow) annotation (Line(points={{5.2,90},
              {0,90},{0,95},{-12,95}}, color={0,0,127}));
      connect(clearance.q_in, concentrationMeasure1.q_in) annotation (Line(
          points={{-22,88},{-26,88},{-26,34},{-94,34},{-94,16}},
          color={107,45,134},
          thickness=1));
      connect(MiddleTorso.solutionVolume, volume.y)
        annotation (Line(points={{72,40},{72,40},{89,40}}, color={0,0,127}));
      connect(UpperTorso.solutionVolume, volume1.y)
        annotation (Line(points={{72,70},{89,70}}, color={0,0,127}));
      connect(LowerTorso.solutionVolume, volume2.y)
        annotation (Line(points={{72,12},{89,12}}, color={0,0,127}));
      connect(Bladder.solutionVolume, volume3.y) annotation (Line(points={{72,
              -10},{72,-10},{89,-10}}, color={0,0,127}));
      connect(plasma.solutionVolume, volume4.y) annotation (Line(points={{-70,
              38},{-70,38},{-70,52},{-69,52}}, color={0,0,127}));
      connect(synthesis.q_out, molarFlowMeasure.q_in) annotation (Line(
          points={{-60,-50},{-54,-50}},
          color={107,45,134},
          thickness=1));
      connect(molarFlowMeasure.q_out, plasma.q_out) annotation (Line(
          points={{-34,-50},{-26,-50},{-26,34},{-66,34}},
          color={107,45,134},
          thickness=1));
      connect(plasma.q_out, molarFlowMeasure1.q_in) annotation (Line(
          points={{-66,34},{-26,34},{-26,-50},{-22,-50}},
          color={107,45,134},
          thickness=1));
      connect(molarFlowMeasure1.q_out, degradation.q_in) annotation (Line(
          points={{-2,-50},{-2,-50},{0,-50}},
          color={107,45,134},
          thickness=1));
      connect(add.u1, molarFlowMeasure1.molarFlowRate) annotation (Line(points=
              {{-1,-62},{-12,-62},{-12,-58}}, color={0,0,127}));
      connect(add.u2, molarFlowMeasure.molarFlowRate) annotation (Line(points={
              {-1,-68},{-1,-68},{-44,-68},{-44,-58}}, color={0,0,127}));
      connect(add.y, proteinDivision1.totalProteins) annotation (Line(points={{
              10.5,-65},{10.5,-65},{20,-65}}, color={0,0,127}));
      connect(proteinDivision1.albumin, proteinCharge.u) annotation (Line(
            points={{30,-62},{40,-62},{40,-61}}, color={0,0,127}));
      connect(pH.y, proteinCharge.pH) annotation (Line(points={{29,-80},{40,-80},
              {40,-79}}, color={0,0,127}));
      connect(proteinCharge.port_a, acidBaseBuffers.port_a) annotation (Line(
          points={{59,-70},{81,-70},{81,-70}},
          color={107,45,134},
          thickness=1));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}})));
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
    equation
    //  ProteinsMassConcentration2Concentration(c.u*1000) = q_in.conc;
      c.u = q_in.conc;
      q_in.q = DegradationBasic * c.val;
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
      Modelica.Blocks.Interfaces.RealInput AlbuminDifferenceMolarFlow
        annotation (Placement(transformation(extent={{-120,70},{-80,110}})));
      constant Real AlbMolarMass( final unit = "g/mol")= 66000;
      FiggeFencl3Detailed figgeFencl3Base(
        pH = pH,
        pCO2=40,
        Pi=1.15,
        alb=u*AlbMolarMass)
        annotation (Placement(transformation(extent={{-20,0},{0,20}})));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}})));
                //  parameter Real alb (unit = "g/dl")= 4.4;
    equation
      port_a.q = figgeFencl3Base.atch;
    end ProteinCharge;

    model AcidBaseBuffers

      Physiolibrary.Chemical.Interfaces.ChemicalPort_a port_a
        annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      Physiolibrary.Types.RealIO.pHOutput pH
        annotation (Placement(transformation(extent={{-80,-100},{-100,-80}})));
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

    end ProteinDivision;
  end AlbuminBorderFlux;
  annotation (uses(Modelica(version="3.2.1"),
      Physiomodel(version="0.2.29"),
      Physiolibrary(version="2.3.1")));
end FullBloodAcidBase;
