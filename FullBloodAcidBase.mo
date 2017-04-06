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
      parameter Real pH0(unit = "mmHg") = 7.4;
      input Real pCO20(unit = "mmHg") = 40;
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

  annotation (uses(
      Physiomodel(version="0.2.29"),
      Physiolibrary(version="2.3.1"),
      Modelica(version="3.2.2")));
end FullBloodAcidBase;
