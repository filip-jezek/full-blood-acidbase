package FullBloodAcidBase

  package PlasmaElectrochemical

    package Auxilliary

      partial model PlasmaElectrochemicalBase
        //   Real SID(displayUnit="meq/l") "Strong ion difference. Normal value 39";
        //   Real pCO2(displayUnit="mmHg") "CO2 partial pressure. Normal value 40";
        //   Real Pi(unit="mmol/l") "Total phosphate. Normal value 1.15";
        //   Real alb(unit="g/dl") "Albumin concentration. Normal value 4.4";

        input Real SID(displayUnit="meq/l")
          "Strong ion difference. Normal value 39";
        input Real pCO2(displayUnit="mmHg")
          "CO2 partial pressure. Normal value 40";
        input Real Pi(unit="mmol/l") "Total phosphate. Normal value 1.15";
        input Real alb(unit="g/dl") "Albumin concentration. Normal value 4.4";
        // Real Alb(unit="mmol/l") = alb*10/66500*1000;
        Real pH(start=10, unit="1");
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end PlasmaElectrochemicalBase;

      model FiggeFencl3
        "Base class for plasma acidbase after Figge and Fencl 3.0, missing inputs for SID, PCO2, total Pi and total albumin. Not meant to be run independently. Implemented after [Figge, 2013]"
        extends
          FullBloodAcidBase.PlasmaElectrochemical.Auxilliary.PlasmaElectrochemicalBase;

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
        constant Real kw=0.000000000000044;

        constant Real Kc1=0.0000000000244
          "Kc1 is derived from the parameters in the Henderson-Hasselbalch equation. pK = 6.1; a = 0.230 mM / kPa; 1 Torr = 0.13332236842105 kPa. The value of Kc1 is 2.44E-11 (Eq / L)^2 / Torr.";
        constant Real Kc2=0.00000000011
          "Kc2 is calculated from Harned and Scholes (1941) for 37 degrees C and ionic strength 0.15 M. The value of Kc2 is 5.5E-11 mol / L x 2 = 1.1E-10 Eq / L.";

        constant Real K1=0.0122
          "K1, K2, and K3 for the phosphoric acid - phosphate system are from Sendroy and Rem: Hastings (1927). pK1 = 1.915; pK2 = 6.66; pK3 = 11.78.";
        constant Real K2=0.000000219
          "K1, K2, and K3 for the phosphoric acid - phosphate system are from Sendroy and Hastings (1927). pK1 = 1.915; pK2 = 6.66; pK3 = 11.78.";
        constant Real K3=0.00000000000166
          "K1, K2, and K3 for the phosphoric acid - phosphate system are from Sendroy and Hastings (1927). pK1 = 1.915; pK2 = 6.66; pK3 = 11.78.";

        //Rem: Enter desired values for SID, PCO2, [ Pi tot ], and [ Albumin ] in the next four lines.
        //   input Real SID(displayUnit="meq/l") "Strong ion difference. Normal value 39";
        //   input Real pCO2(displayUnit="mmHg") "CO2 partial pressure. Normal value 40";
        //   input Real Pi(unit="mmol/l") "Total phosphate. Normal value 1.15";
        //   input Real alb(unit="g/dl") "Albumin concentration. Normal value 4.4";

      public
        Real H(displayUnit="eq/l") = 10^(-pH);
        Real HO(displayUnit="eq/l") = kw/H;
        Real HCO3(displayUnit="eq/l") = Kc1*pCO2/H;
        Real CO3(displayUnit="eq/l") = Kc2*HCO3/H;

      protected
        Real FNX=K1*H^2 + 2*K1*K2*H + 3*K1*K2*K3;
        Real FNY=H^3 + K1*H^2 + K1*K2*H + K1*K2*K3;
        Real FNZ=FNX/FNY;
      public
        Real P(displayUnit="meq/l") = Pi*FNZ;
        Real Netcharge=SID + 1000*(H - HO - HCO3 - CO3) - P;

        Real NB=0.4*(1 - (1/(1 + 10^(pH - 6.9))))
          "NB accounts for histidine pK shift due to the NB transition";

        // constant Real albuminResidues[:] = cat(1,{-1 /*cysteine */,-98/*glutamic acid*/,-18/*tyrosine*/,+24/*arginine */, /* lysine >>>*/ 2, 2, 2, 2, 1, 50} ,ones(16) /*histidine residues*/,/* amino terminus and carboxyl terminus*/{1, 1});
      protected
        Real albConversion=1000*10*alb/66500;
        constant Real albuminResidues[:]=cat(
                  1,
                  {-1,-98,-18,+24,2,2,2,2,1,50},
                  ones(16),
                  {1,-1});
        // Real albuminPks[:] = {8.5 /* CYST*/,3.9 /* GLUT*/,11.7 /* TYR*/,12.5 /* ARG*/,/*LYS >>>*/5.8, 6.15, 7.51, 7.685, 7.86, 10.3,/*HIST>>>*/7.12 - NB, 7.22 - NB, 7.1 - NB, 7.49 - NB, 7.01 - NB, 7.31, 6.75, 6.36, 4.85, 5.76, 6.17, 6.73, 5.82, 5.1, 6.7, 6.2, 8/* amino terminus */,3.1 /*carboxyl terminus*/};
        Real albuminPks[:]={8.5,3.9,11.7,12.5,5.8,6.15,7.51,7.685,7.86,10.3,
            7.12 - NB,7.22 - NB,7.1 - NB,7.49 - NB,7.01 - NB,7.31,6.75,6.36,
            4.85,5.76,6.17,6.73,5.82,5.1,6.7,6.2,8,3.1};
        Real albChrg[n](each displayUnit="meq/l") "charge of albumin per unit";
        constant Integer n=size(albuminResidues, 1);
      public
        Real atch=sum(albChrg)*albConversion
          "albumin total charge. Normal value -12.2678";

      equation
        Netcharge + atch = 0;

        for i in 1:n loop
          if albuminResidues[i] > 0 then
            // positive charge
            albChrg[i] = albuminResidues[i]*(1/(1 + 10^(pH - albuminPks[i])));
          else
            // negative charge
            albChrg[i] = albuminResidues[i]*(1/(1 + 10^(-pH + albuminPks[i])));
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
        annotation (Documentation(info="<html>
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

      model FiggeFencl3Detailed
        "Extension for investigation of detailed albumin balance"
        extends Auxilliary.FiggeFencl3;
      protected
        Real albPositivePart[n], albNegativePart[n];
        Real albTotalPlusPart[n], albTotalMinusPart[n];
      public
        Real HCO3mEqL=HCO3*1000;
        Real albHAPlus=albConversion*sum(albPositivePart) "A0 + H+ = HA+";
        Real albAMinus=-albConversion*sum(albNegativePart) "A- + H+ = HA0";
        Real albA0=(ATotPlus - albHAPlus) "A0 + H+ = HA+";
        Real albHA0=(ATotMinus - albAMinus) "HA0 = A- + H+";
        Real ATotPlus=albConversion*sum(albTotalPlusPart)
          "Part of albumin, which could be positive";
        Real ATotMinus=-albConversion*sum(albTotalMinusPart)
          "Part of albumin, which could be negative.";

        Real ach=sum(albChrg);
        Real atch0=-12.2678 "to demonstrate the low buffer strength of albumin";
        Real atot1=albHAPlus + albA0;
        Real atot2=albAMinus + albHA0;
        Real test=atot1 + atot2;

        Real barGraphAlb1=albHA0;
        Real barGraphAlb2=albHA0 + albAMinus;
        Real barGraphAlb3=albHA0 + albAMinus + albHAPlus;
        Real barGraphAlb4=albHA0 + albAMinus + albHAPlus + albA0;
        Real barGraphHCO3Alb1=P;
        Real barGraphHCO3Alb2=P - atch;
        Real barGraphHCO3Alb3=P - atch + 1000*HCO3;
        discrete Real dHCO3s, dHCO3e, dHCO3;
        discrete Real dAlbs, dAlbe, dAlb;
        Real diff=dHCO3 - dAlb;
        Real outofcobntrol=albHAPlus - albAMinus;
        parameter Real p=1;
        parameter Real a=1;
        parameter Real b=1;
        parameter Real c=1;
        Real AlbXMinus=-(0.148*pH - 0.818)*(alb*10) "Total albumin charge";
        Real PXminus=-(0.309*pH - 0.469)*Pi;
        Real Pminus=-P;
        Real totalDiff=AlbXMinus - atch + PXminus - Pminus;

      equation

        // POSSIBLE OPENMODELICA INCOMPATIBILITY
        // IN OPENMODELICA USE ELSEWHEN instead of end when; and when terminal.. lines
        when initial() then
          dHCO3s = 1000*HCO3;
          dAlbs = albHAPlus - albAMinus - p*P - 1000*CO3*a - 1000*HO*b + 1000*H
            *c;
        end when;
        when terminal() then
          // elsewhen terminal() then
          dHCO3e = 1000*HCO3;
          dHCO3 = dHCO3e - dHCO3s;
          dAlbe = albHAPlus - albAMinus - p*P - 1000*CO3*a - 1000*HO*b + 1000*H
            *c;
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
        annotation (Documentation(info="<html>
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

      model Fencl
        "Simple albumin and Pi calculations, according to [Fencl, 2000]"
        extends
          FullBloodAcidBase.PlasmaElectrochemical.Auxilliary.PlasmaElectrochemicalBase;

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
        constant Real kw=0.000000000000044;

        constant Real Kc1=0.0000000000244
          "Kc1 is derived from the parameters in the Henderson-Hasselbalch equation. pK = 6.1; a = 0.230 mM / kPa; 1 Torr = 0.13332236842105 kPa. The value of Kc1 is 2.44E-11 (Eq / L)^2 / Torr.";
        constant Real Kc2=0.00000000011
          "Kc2 is calculated from Harned and Scholes (1941) for 37 degrees C and ionic strength 0.15 M. The value of Kc2 is 5.5E-11 mol / L x 2 = 1.1E-10 Eq / L.";

        constant Real K1=0.0122
          "K1, K2, and K3 for the phosphoric acid - phosphate system are from Sendroy and Rem: Hastings (1927). pK1 = 1.915; pK2 = 6.66; pK3 = 11.78.";
        constant Real K2=0.000000219
          "K1, K2, and K3 for the phosphoric acid - phosphate system are from Sendroy and Hastings (1927). pK1 = 1.915; pK2 = 6.66; pK3 = 11.78.";
        constant Real K3=0.00000000000166
          "K1, K2, and K3 for the phosphoric acid - phosphate system are from Sendroy and Hastings (1927). pK1 = 1.915; pK2 = 6.66; pK3 = 11.78.";

        // Real H(displayUnit="eq/l") = 10^(-pH);
      public
        Real HCO3(displayUnit="mmol/l") = 1000*Kc1*pCO2/10^(-pH);
        Real P(displayUnit="meq/l") = -Pi*(0.309*pH - 0.469);
        Real atch=-(alb*10)*(0.123*pH - 0.631) "albumin total charge";
      equation
        SID + P + atch - HCO3 = 0;
        annotation (Documentation(info="<html>
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
      end Fencl;

      model Wolf
        "Simple albumin and Pi calculations, according to [Wolf, 2011]"
        extends
          FullBloodAcidBase.PlasmaElectrochemical.Auxilliary.PlasmaElectrochemicalBase;
        import Modelica.SIunits.*;

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

        constant Real Vp0ByVp=0.98468;
        constant Real fpw=0.9406 "fraction plasma water - Vpw/Vp";


        parameter Concentration Na=140;
        parameter Concentration K=4;
        parameter Concentration Ca=2.4;
        parameter Concentration Mg=0.8;
        // already inherited
        //  parameter Concentration Pi=1.2;
        Real AlbPwGPerL=alb*10*Vp0ByVp "g/lpw";
        Concentration Cl;
        Concentration Alb=AlbPwGPerL/66.5 "mmol/Lpw";
        parameter Concentration im=11.87;
        parameter Concentration Lac=1.5;
        Concentration SO4pw=0.33/fpw "mmol/Lpw";
        Concentration volume_c[cont]={Na,K,Ca,Mg,Cl,Pi,Alb/Vp0ByVp,im,Lac}
          "concentration in one liter";
        Concentration water_c[cont]=volume_c/fpw*Vp0ByVp
          "Actual concentration recalculated to water fraction in one liter (initially 0.96)";

        Real CaPPBound=(3.98)/(3.98 + Hh)*water_c[cont.Ca] "1.257";
        Real MgPPBound=CaPPBound/2;
        Real ZCaBindPerAlb=CaPPBound/water_c[cont.Alb];
        Real ZClBindPerAlb=6.46/(6.46 + Hh)*6 + 4 "Anstey";
        Real CaIon=water_c[cont.Ca]*2 - CaPPBound "mEq/Lpw";
        Real MgIon=water_c[cont.Mg]*2 - MgPPBound "mEq/Lpw";
        //charge on inpermeable solutes
        Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87));
        Real ZFigge=(-10.65) - 16*(10^(pH - 7.418)/(1 + 10^(pH - 7.418)));
        Real ZAlbBnd=-ZClBindPerAlb + ZCaBindPerAlb + ZCaBindPerAlb/2;
        Real ZAlb=ZFigge + ZAlbBnd;
        Real Hh=H/fpw*1e8;

        //   parameter Real Zim ;//= -5.3 "Charge of ALL impermeable solutes";//test
        // Real Zim = -5.3 "Charge of ALL impermeable solutes";//test
        parameter Real Zim=-0.6349 "Charge of ALL impermeable solutes";

        //   Real fi[cont]={0.93,0.93,1,1,0.93,0.93,1,1,-9999};
        // Alb has osmotic coefficient of 1???
        //   Real CaOsm=water_c[cont.Ca] - CaPPBound*0.5;
        //   Real MgOsm=water_c[cont.Mg] - MgPPBound*0.5;
        //   Real OsmPart=SO4pw + CaOsm + MgOsm + sum(water_c[{cont.Na,cont.K,cont.Cl,cont.Lac}])
        //        + water_c[cont.Pi] + HCO3 + CO3;
        //   Real Osm=OsmPart*0.93 + permeableParticles " = 284.084";

        //   Real Z[cont]={1,1,2,2,-1,ZPi,ZAlb,0,1};
        //   Real Chrgs[cont]=water_c .* Z;

        Real atch=water_c[cont.Alb]*ZFigge;
        Real charge=Zim + SID + water_c[cont.Alb]*ZAlb + water_c[cont.Pi]*ZPi
             - HCO3 - 2*CO3 - SO4pw*2;

        Concentration HCO3=0.0306*pCO2mmHg*10^(pH - 6.103)/fpw;
        Concentration CO3=HCO3*10^(pH - 10.2);
        //  Real pH(start=7.37); = -log10(H);


        Real pCO2mmHg(unit="1") = pCO2;
        Real H(
          start=10^(-7.37),
          min=0,
          max=1) = 10^(-pH);
        Real Hw=H/fpw;
        Concentration permeableParticles=(5 + 5)/fpw
          "glucose and urea concentration in PLasma water";
      equation
        //  SID = HCO3 + P + atch;
        charge = 0;
        SID = sum(water_c[{cont.Na,cont.K}]) - water_c[cont.Cl] + 2*sum(water_c[
          {cont.Ca,cont.Mg}]) - CaPPBound - MgPPBound - water_c[cont.Lac];
        annotation (Documentation(info="<html>
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
      end Wolf;

      model FiggeFenclNSID "Calculation of normal SID using plasma figge fencl"
        parameter Real pH0(displayUnit="mmHg") = 7.4;
        input Real pCO20(displayUnit="mmHg") = 40;
        input Real Pi0(unit="mmol/l");
        //= 1.15;
        input Real alb0(unit="g/dl");
        //= 4.4;
        Real SID;
        // = figgeFencl3.SID;

        Auxilliary.FiggeFencl3 figgeFencl3(
          pH=pH0,
          pCO2=pCO20,
          Pi=Pi0,
          alb=alb0,
          SID=SID)
          annotation (Placement(transformation(extent={{-58,0},{-38,20}})));
      end FiggeFenclNSID;
    end Auxilliary;

    model PlasmaFiggeFencl3
      "Implementation of FF3 [Figge, 2013] together with NSID"
      extends FullBloodAcidBase.PlasmaElectrochemical.Auxilliary.FiggeFencl3;

      constant Real pH0(displayUnit="mmHg") = 7.4;
      constant Real pCO20(displayUnit="mmHg") = 40;
      Real NSID;
      //  input Real BE=SID - NSID;
      //input Real Pi;
      //input Real alb;
    protected
      FullBloodAcidBase.PlasmaElectrochemical.Auxilliary.FiggeFencl3
        normalPlasmaSID(
        pCO2=pCO20,
        Pi=Pi,
        alb=alb,
        SID=NSID)
        annotation (Placement(transformation(extent={{-58,0},{-38,20}})));


    equation

      normalPlasmaSID.pH = pH0;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end PlasmaFiggeFencl3;

    model PlasmaFencl
      "Implementation of Fencl [Fencl, 2000] together with NSID"
      extends FullBloodAcidBase.PlasmaElectrochemical.Auxilliary.Fencl;

      constant Real pH0(displayUnit="mmHg") = 7.4;
      constant Real pCO20(displayUnit="mmHg") = 40;
      Real NSID;
      //  Real BE=SID - NSID;
      //   input Real Pi;
      //   input Real alb;
    protected
      FullBloodAcidBase.PlasmaElectrochemical.Auxilliary.Fencl normalPlasmaSID(
        pCO2=pCO20,
        Pi=Pi,
        alb=alb,
        SID=NSID)
        annotation (Placement(transformation(extent={{-58,0},{-38,20}})));

    equation

      normalPlasmaSID.pH = pH0;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end PlasmaFencl;

    model PlasmaWolf
      "Implementation of Wolf model of plasma [Wolf, 2011] together with NSID"
      extends FullBloodAcidBase.PlasmaElectrochemical.Auxilliary.Wolf;

      constant Real pH0(displayUnit="mmHg") = 7.4;
      constant Real pCO20(displayUnit="mmHg") = 40;
      Real NSID;
      //  Real BE=SID - NSID;
      //   input Real Pi;
      //   input Real alb;
    protected
      FullBloodAcidBase.PlasmaElectrochemical.Auxilliary.Wolf normalPlasmaSID(
        pCO2=pCO20,
        Pi=Pi,
        alb=alb,
        SID=NSID)
        annotation (Placement(transformation(extent={{-58,0},{-38,20}})));

    equation

      normalPlasmaSID.pH = pH0;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end PlasmaWolf;

    package Tests
      extends Modelica.Icons.ExamplesPackage;
      model Plasma_pCO2 "Change of pH while changing pCO2"
        extends Modelica.Icons.Example;

        FullBloodAcidBase.PlasmaElectrochemical.PlasmaFiggeFencl3 ff(
          SID=ff.NSID + BE,
          pCO2=pCO2,
          Pi=Pi,
          alb=alb)
          annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
        FullBloodAcidBase.PlasmaElectrochemical.PlasmaFencl fe(
          SID=fe.NSID + BE,
          pCO2=pCO2,
          Pi=Pi,
          alb=alb)
          annotation (Placement(transformation(extent={{-60,20},{-40,40}})));

        FullBloodAcidBase.PlasmaElectrochemical.PlasmaWolf wo(
          SID=wo.NSID + BE,
          pCO2=pCO2,
          Pi=Pi,
          alb=alb)
          annotation (Placement(transformation(extent={{-60,20},{-40,40}})));

        //parameter Real SID=39;
        //   Real SID=39 - 15 + time*30;
        //   Real pCO2=40;
        parameter Real BE=0;
        Real pCO2=time*40 + 20;
        parameter Real Pi=1.15;
        parameter Real alb=4.4;
        output Real pHFF=ff.pH;
        output Real pHFe=fe.pH;
        output Real pHWo=wo.pH;
      equation
        //  BE =
        annotation (experiment(Tolerance=1e-006),
            __Dymola_experimentSetupOutput);
      end Plasma_pCO2;

      model Plasma_BE "Change of pH while changing pCO2"
        extends Modelica.Icons.Example;

        FullBloodAcidBase.PlasmaElectrochemical.PlasmaFiggeFencl3 ff(
          SID=ff.NSID + BE,
          pCO2=pCO2,
          Pi=Pi,
          alb=alb)
          annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
        FullBloodAcidBase.PlasmaElectrochemical.PlasmaFencl fe(
          SID=fe.NSID + BE,
          pCO2=pCO2,
          Pi=Pi,
          alb=alb)
          annotation (Placement(transformation(extent={{-60,20},{-40,40}})));

        FullBloodAcidBase.PlasmaElectrochemical.PlasmaWolf wo(
          SID=wo.NSID + BE,
          pCO2=pCO2,
          Pi=Pi,
          alb=alb)
          annotation (Placement(transformation(extent={{-60,20},{-40,40}})));

        //parameter Real SID=39;
        //   Real SID=39 - 15 + time*30;
        //   Real pCO2=40;
        Real BE=-15 + time*30;
        parameter Real pCO2=40;
        parameter Real Pi=1.15;
        parameter Real alb=4.4;
        output Real pHFF=ff.pH;
        output Real pHFe=fe.pH;
        output Real pHWo=wo.pH;
      equation
        //  BE =
        annotation (experiment(Tolerance=1e-006),
            __Dymola_experimentSetupOutput);
      end Plasma_BE;
    end Tests;
  end PlasmaElectrochemical;

  package FullBloodEmpirical
    package Auxiliary
      partial model FullBloodEmpiricalBase
        Real Hb(unit="g/dl") = Hct*33.34;
        input Real Hct(unit="1") = 15/33.34;
        input Real BEox(displayUnit="meq/l") = 0;
        input Real sO2(unit="1") = 1;
        input Real pCO2(displayUnit="mmHg") = 40;
        Real pH(start=7.4);
      end FullBloodEmpiricalBase;
    end Auxiliary;

    model SAVanSlyke77
      "original Van Slyke equation according to SA 1977 DOI: 10.3109/00365517709098927"
      extends Auxiliary.FullBloodEmpiricalBase;
      Real a(unit="mmol/l") "bicarbonate concentration in plasma/(mmol/l)";
      Real b(unit="mmol/l") = Hb*0.6206
        "hemoglobin concentration in blood/(mmol/l)";
      Real c(unit="1") = pH "pH of plasma at 37 degrees C";
      Real d(displayUnit="mEq/l") = BEox - 0.3*(1 - sO2)
        "base excess concentration in blood/(mmol/l)";

    equation
      a - 24.4 = -(2.3*b + 7.7)*(c - 7.40) + d/(1 - 0.023*b);
      pH = 6.1 + log10(a/(0.230*pCO2*133/1000));
    end SAVanSlyke77;

    model Zander1995 "Implementation according to [Zander, 1995]"
      extends Auxiliary.FullBloodEmpiricalBase;
      Real BE=(1 - 0.0143*Hb)*((0.0304*pCO2*(10^(pH - 6.1)) - 24.26) + (9.5 +
          1.63*Hb)*(pH - 7.4)) - 0.2*Hb*(1 - sO2);
    equation
      BE = BEox + 0.2*(1 - sO2)*Hb;
    end Zander1995;

    model SAVanSlyke
      "The Van Slyke equation by Sigaard-Andersen, taken from http://www.siggaard-andersen.dk/OsaTextbook.htm"
      extends Auxiliary.FullBloodEmpiricalBase;
      //    parameter Real Hct = 15/33.34;
      //    Real BEox = 0;
      //    parameter Real sO2 = 1;
      //    parameter Real pCO2 = 40;

      Real ctH=-(1 - (1 - rc)*fiEB)*((cHCO3 - cHCO30) + buf*(pH - 7.4))
        "concentration of titratable hydrogen ion";

      // Real rc = 0.57; // cHCO3e/cHCO3p from web pages
      Real rc=0.57 - 0.28*dpH - 0.082*dpH^2;
      // said to be from siggaard-andersen, 1974: The acid-base status of the blood
      Real dpH=-(pH - 7.4);

      //parameter Real fiEB = 0.7;//15/33.34;// Hct;// ctHbB/cHBE
      Real fiEB=Hct;
      // ctHbB/cHBE
      Real ctHbE(unit="mmol/l") = 21;
      Real cHCO30(unit="mmol/l") = 24.5;
      Real buf=2.3*ctHb + BetaP;
      Real ctHb=Hb*0.6206
        "Hb*4 SA is considering the monomer of hemoglobin tetramer";
      //Hct*33.34; // recompute the hematocrit to mmol/l - BUT IS THAT COMPATIBLE? TODO!
      Real BetaP=5.8 + 8.0*(cAlb - 0.66);
      // Real BetaP = 7.7 + 8.0 * (cAlb - 0.66);
      Real cAlb(unit="mmol/l") = Alb/6.6;
      //0.66;
      input Real Alb(unit="g/dl") = 4.4;

      Real BE=-ctH;
      Real cHCO3(start=24.5);
    equation
      BEox = BE + 0.3*(1 - sO2);

      pH = 6.1 + log10(cHCO3/(0.230*pCO2*133/1000));

    end SAVanSlyke;

    model Kofr2009 "[Kofranek, 2009]"
      extends Auxiliary.FullBloodEmpiricalBase;

      Real a[:]={996.35 - 10.35*Hb,35.16875 + 0.25875*Hb,-82.41 + 2.01*Hb,-5.27625
           - (5.025e-2)*Hb,121 - Hb,2.625 + 0.025*Hb,-2.556 - 0.0944*Hb,
          13.87634 + 0.186653*Hb + 0.00534936*Hb^2,0.548 - 0.0274*Hb,0.274 -
          0.0137*Hb};

    public
      Real BE=BEox + 0.2*(1 - sO2)*Hb;
      Real Y=a[7] + sqrt(a[8] + a[9]*BE);

    equation
      pH = (a[1]*a[10] + a[2]*Y + (a[3]*a[10] + a[4]*Y)*log10(pCO2))/(a[5]*a[10]
         + a[6]*Y);

    end Kofr2009;

    model SAoriginal
      "Our complex polynomial fit to original SA data, used as a reference"
      extends Auxiliary.FullBloodEmpiricalBase;

    protected
      constant Real pco2BBCoef[:]={2.1125e-009,-640.9926e-009,72.7649e-006,-3.2862e-003,
          -38.1749e-003,8.2352e+000,-97.0551e+000};
      constant Real pco2BECoef[:]={8.3975e-009,-513.9503e-009,3.8105e-006,
          231.6888e-006,-46.5581e-003,353.7105e-003,39.9871e+000};
      constant Real pHBBCoef[:]={40.8936e-012,-13.0063e-009,1.6780e-006,-111.7919e-006,
          4.0776e-003,-67.8274e-003,7.2888e+000};
      constant Real pHBECoef[:]={131.3315e-009,2.5027e-006,175.6144e-006,
          11.9273e-003,7.4001e+000};
    public
      Real pCO2BB(start=96) = pco2BBCoef[1]*BB^6 + pco2BBCoef[2]*BB^5 +
        pco2BBCoef[3]*BB^4 + pco2BBCoef[4]*BB^3 + pco2BBCoef[5]*BB^2 +
        pco2BBCoef[6]*BB + pco2BBCoef[7];
      Real pCO2BE(start=40) = pco2BECoef[1]*BE^6 + pco2BECoef[2]*BE^5 +
        pco2BECoef[3]*BE^4 + pco2BECoef[4]*BE^3 + pco2BECoef[5]*BE^2 +
        pco2BECoef[6]*BE + pco2BECoef[7];
      Real pHBB(start=7) = pHBBCoef[1]*BB^6 + pHBBCoef[2]*BB^5 + pHBBCoef[3]*BB
        ^4 + pHBBCoef[4]*BB^3 + pHBBCoef[5]*BB^2 + pHBBCoef[6]*BB + pHBBCoef[7];
      Real pHBE(start=7) = pHBECoef[1]*BE^4 + pHBECoef[2]*BE^3 + pHBECoef[3]*BE
        ^2 + pHBECoef[4]*BE + pHBECoef[5];

      Real BB=BE + 0.42*Hb + 41.7;
      Real BE=BEox + 0.2*(1 - sO2)*Hb;

    equation
      pH = (log10(pCO2) - log10(pCO2BB))*(pHBB - pHBE)/(log10(pCO2BB) - log10(
        pCO2BE)) + pHBB;

    end SAoriginal;

    package Tests
      extends Modelica.Icons.ExamplesPackage;
      model Blood_pCO2
        "pH dependent on varying pCO2 by different approximations. We take our implementation of SA nomogram as reference. Compare the object's pH"
        extends Modelica.Icons.Example;
        parameter Real Hct=15/33.34;
        parameter Real BEox=0;
        parameter Real sO2=1;
        output Real pHZander=zander1995.pH;
        output Real pHNomogram=sAoriginal.pH;
        output Real pHVanSlyke=sAVanSlyke.pH;
        output Real pHVanSlyke77=sAVanSlyke77.pH;
        input Real pCO2=time*40 + 20;
        FullBloodEmpirical.Zander1995 zander1995(
          pCO2=pCO2,
          BEox=BEox,
          Hct=Hct,
          sO2=sO2)
          annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
        FullBloodEmpirical.Kofr2009 kofr2009(
          pCO2=pCO2,
          BEox=BEox,
          Hct=Hct,
          sO2=sO2)
          annotation (Placement(transformation(extent={{-20,40},{0,60}})));
        FullBloodEmpirical.SAoriginal sAoriginal(
          pCO2=pCO2,
          BEox=BEox,
          Hct=Hct,
          sO2=sO2)
          annotation (Placement(transformation(extent={{34,40},{54,60}})));
        FullBloodEmpirical.SAVanSlyke sAVanSlyke(
          pCO2=pCO2,
          BEox=BEox,
          Hct=Hct,
          sO2=sO2)
          annotation (Placement(transformation(extent={{80,40},{100,60}})));
        FullBloodEmpirical.SAVanSlyke77 sAVanSlyke77(
          pCO2=pCO2,
          BEox=BEox,
          Hct=Hct,
          sO2=sO2)
          annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
        annotation (
          experiment(Tolerance=0.001),
          __Dymola_experimentSetupOutput,
          __Dymola_Commands(file="def.mos" "def"));
      end Blood_pCO2;

      model Blood_BE
        "pH dependent on varying pCO2 by different approximations. We take our implementation of SA nomogram as reference. Compare the object's pH"
        extends Modelica.Icons.Example;
        parameter Real Hct=15/33.34;
        Real BEox=-15 + time*30;
        parameter Real sO2=1;
        output Real pHZander=zander1995.pH;
        output Real pHNomogram=sAoriginal.pH;
        output Real pHVanSlyke=sAVanSlyke.pH;
        output Real pHVanSlyke77=sAVanSlyke77.pH;
        parameter Real pCO2=40;
        //time*40 + 20;
        FullBloodEmpirical.Zander1995 zander1995(
          pCO2=pCO2,
          BEox=BEox,
          Hct=Hct,
          sO2=sO2)
          annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
        FullBloodEmpirical.Kofr2009 kofr2009(
          pCO2=pCO2,
          BEox=BEox,
          Hct=Hct,
          sO2=sO2)
          annotation (Placement(transformation(extent={{-20,40},{0,60}})));
        FullBloodEmpirical.SAoriginal sAoriginal(
          pCO2=pCO2,
          BEox=BEox,
          Hct=Hct,
          sO2=sO2)
          annotation (Placement(transformation(extent={{34,40},{54,60}})));
        FullBloodEmpirical.SAVanSlyke sAVanSlyke(
          pCO2=pCO2,
          BEox=BEox,
          Hct=Hct,
          sO2=sO2)
          annotation (Placement(transformation(extent={{80,40},{100,60}})));
        FullBloodEmpirical.SAVanSlyke77 sAVanSlyke77(
          pCO2=pCO2,
          BEox=BEox,
          Hct=Hct,
          sO2=sO2)
          annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
        annotation (
          experiment(Tolerance=0.001),
          __Dymola_experimentSetupOutput,
          __Dymola_Commands(file="def.mos" "def"));
      end Blood_BE;
    end Tests;
    annotation (Documentation(info="<html>
<h4>References:</h4>
<p>Fencl, V., Jabor, A., Kazda, A., &amp; Figge, J. (2000). Appendix of Diagnosis of Metabolic Acid&ndash;Base Disturbances in Critically Ill Patients. American Journal of Respiratory and Critical Care Medicine, 162(6), 2246&ndash;2251.</p>
<p>Figge, J. (27 October, 2013). The Figge-Fencl Quantitative Physicochemical Model of Human Acid-Base Physiology (Version 3.0). Retrieved June 22, 2016, from http://www.figge-fencl.org/model.html</p>
<p>Kofr&aacute;nek, J. (2009). KOMPLEXN&Iacute; MODEL ACIDOBAZICK&Eacute; ROVNOV&Aacute;HY KRVE. Medsoft 2009. Retrieved from http://www.physiome.cz/references/Medsoft09a.pdf</p>
<p>Raftos, J. E., Bulliman, B. T., &amp; Kuchel, P. W. (1990). Evaluation of an electrochemical model of erythrocyte pH buffering using 31P nuclear magnetic resonance data. The Journal of General Physiology, 95(6), 1183&ndash;1204.</p>
<p>Siggaard-Andersen, O. (1977). The van Slyke equation. Scandinavian Journal of Clinical and Laboratory Investigation. Supplementum, 146, 15&ndash;20.</p>
<p>Siggaard-Andersen, O., &amp; Fogh-Andersen, N. (1995). Base excess or buffer base (strong ion difference) as measure of a non-respiratory acid-base disturbance. Acta Anaesthesiologica Scandinavica. Supplementum, 107, 123&ndash;128.</p>
<p>Wolf, M. B. (2013). Whole body acid-base and fluid-electrolyte balance: a mathematical model. American Journal of Physiology. Renal Physiology, 305(8), F1118&ndash;31.</p>
<p>Wolf, M. B. (2015). Comprehensive diagnosis of whole-body acid-base and fluid-electrolyte disorders using a mathematical model and whole-body base excess. Journal of Clinical Monitoring and Computing, 29(4), 475&ndash;490.</p>
<p>Wolf, M. B., &amp; DeLand, E. C. (2011). A comprehensive, computer-model-based approach for diagnosis and treatment of complex acid&ndash;base disorders in critically-ill patients. Journal of Clinical Monitoring and Computing, 25(6), 353&ndash;364.</p>
<p>Wolf, M. B., &amp; Deland, E. C. (2011). A mathematical model of blood-interstitial acid-base balance: application to dilution acidosis and acid-base status. Journal of Applied Physiology, 110(4), 988&ndash;1002.</p>
<p>Zander, R. (1995). Die korrekte Bestimmung des Base Excess (BE, mmol/l) im Blut. AINS-An&auml;sthesiologie- Intensivmedizin- Notfallmedizin- Schmerztherapie, 30(S 1), S36&ndash;S38.</p>
</html>"));
  end FullBloodEmpirical;

  package FullBloodCombined

    model CombinedModel
      "Test combined model of Figge-fencl plasma and SA full hemoatocrite"
      replaceable PlasmaElectrochemical.PlasmaFiggeFencl3 plasma(
        SID=SID,
        pCO2=pCO2,
        Pi=Pi,
        alb=alb) constrainedby
        FullBloodAcidBase.PlasmaElectrochemical.Auxilliary.PlasmaElectrochemicalBase(
        SID=SID,
        pCO2=pCO2,
        Pi=Pi,
        alb=alb)
        annotation (Placement(transformation(extent={{-20,20},{0,40}})));

      replaceable FullBloodEmpirical.SAoriginal fullErythrocyte(
        Hct=1,
        BEox=BEe,
        pCO2=pCO2) constrainedby
        FullBloodAcidBase.FullBloodEmpirical.Auxiliary.FullBloodEmpiricalBase(
        Hct=1,
        BEox=BEe,
        pCO2=pCO2)
        annotation (Placement(transformation(extent={{20,20},{40,40}})));
      constant Real fullHb=33.34;
      input Real Hb=15;
      input Real Hct=Hb/fullHb;

      Real BEp(unit="meq/l") = BE - mHCO3/(1 - Hct);
      Real BEe(unit="meq/l") = BE + mHCO3/Hct;

      Real mHCO3;

      Real SID;
      input Real BE=0;
      input Real pCO2=40;
      input Real Pi=1.15;
      input Real alb=4.4;
      output Real pH=fullErythrocyte.pH;
    equation
      plasma.pH = fullErythrocyte.pH;
      BEp = SID - plasma.NSID;
      annotation (experiment(
          StopTime=1,
          __Dymola_NumberOfIntervals=500,
          Tolerance=1e-003));
    end CombinedModel;

    package comparisson

      package Auxiliary

        model ResultSetAtBE
          parameter Real Pi=1.15;
          parameter Real alb=4.4;
          parameter Real BE=0;
          parameter Real Hb=15;
          output Real pHSABlood=normalBloodSA.pH;
          output Real pHComboBlood=Normal.full_blood_combo.pH;
          output Real pHFFplasma=Normal.FF_plasma_only.pH;
          output Real pHSAplasma=plasmaSA.pH;
          Real logpCO2=log10(pCO2);
          input Real pCO2=time*40 + 20;

          FullBloodAcidBase.FullBloodCombined.comparisson.Auxiliary.SetAtAlb
            LowAlb(
            BE=BE,
            pCO2=pCO2,
            alb=alb/2,
            Hb=Hb,
            Pi=Pi) "Low Alb"
            annotation (Placement(transformation(extent={{20,20},{40,40}})));

          FullBloodAcidBase.FullBloodCombined.comparisson.Auxiliary.SetAtAlb
            Normal(
            BE=BE,
            pCO2=pCO2,
            alb=alb,
            Hb=Hb,
            Pi=Pi) "Normal model set"
            annotation (Placement(transformation(extent={{20,20},{40,40}})));

          FullBloodAcidBase.FullBloodCombined.comparisson.Auxiliary.SetAtAlb
            HighAlb(
            BE=BE,
            pCO2=pCO2,
            alb=alb*2,
            Hb=Hb,
            Pi=Pi) "High Alb"
            annotation (Placement(transformation(extent={{20,20},{40,40}})));
          FullBloodEmpirical.SAoriginal plasmaSA(
            Hb=0,
            BEox=BE,
            pCO2=pCO2) "SA original with Hct 0 for comparison with FF"
            annotation (Placement(transformation(extent={{60,20},{80,40}})));

          FullBloodEmpirical.SAoriginal normalBloodSA(
            Hb=Hb,
            BEox=BE,
            pCO2=pCO2) "SA original for comparison only"
            annotation (Placement(transformation(extent={{60,20},{80,40}})));

        end ResultSetAtBE;

        model SetAtAlb

          CombinedModel full_blood_combo(
            BE=BE,
            pCO2=pCO2,
            alb=alb,
            Hb=Hb,
            Pi=Pi)
            annotation (Placement(transformation(extent={{20,20},{40,40}})));
          FullBloodAcidBase.PlasmaElectrochemical.PlasmaFiggeFencl3
            FF_plasma_only(
            SID=FF_plasma_only.NSID + BE,
            pCO2=pCO2,
            Pi=Pi,
            alb=alb) "Plasma model acc to FiggeFencl for comparison only"
            annotation (Placement(transformation(extent={{60,-20},{80,0}})));

          FullBloodEmpirical.SAoriginal normalBloodSA(
            Hb=Hb,
            BEox=BE,
            pCO2=pCO2) "SA original for comparison only"
            annotation (Placement(transformation(extent={{60,20},{80,40}})));
          FullBloodAcidBase.FullBloodEmpirical.SAVanSlyke sAVanSlyke(
            pCO2=pCO2,
            BEox=BE,
            Hct=Hb/33.34,
            sO2=1,
            Alb=alb)
            annotation (Placement(transformation(extent={{-14,40},{6,60}})));
          Wolf_full_blood wolf_full_blood(
            pCO2=pCO2,
            BE=BE,
            AlbP=alb/AlbMW*1000*10/0.94);
          input Real BE, pCO2;
          input Real alb, Pi;
          input Real Hb;
        protected
          constant Real AlbMW=66463;
          //   PlasmaElectrochemical.O.FiggeFenclNSID normalPlasma(
          //     pH0=7.4,
          //     pCO20=40,
          //     alb0=alb,
          //     Pi0=Pi) annotation (Placement(transformation(extent={{-80,40},{-60,60}})));

        end SetAtAlb;

        model CombinedModelBase
          "Test combined model of generalized physicochemical plasma and SA full hematocrite"
          replaceable PlasmaElectrochemical.PlasmaFiggeFencl3 plasma(
            SID=SID,
            pCO2=pCO2,
            Pi=Pi,
            alb=alb) constrainedby
            FullBloodAcidBase.PlasmaElectrochemical.Auxilliary.PlasmaElectrochemicalBase(
            SID=SID,
            pCO2=pCO2,
            Pi=Pi,
            alb=alb)
            annotation (Placement(transformation(extent={{-20,20},{0,40}})));

          replaceable FullBloodEmpirical.SAoriginal fullErythrocyte(
            Hct=1,
            BEox=BEe,
            pCO2=pCO2) constrainedby
            FullBloodAcidBase.FullBloodEmpirical.Auxiliary.FullBloodEmpiricalBase(
            Hct=1,
            BEox=BEe,
            pCO2=pCO2)
            annotation (Placement(transformation(extent={{20,20},{40,40}})));
          constant Real fullHb=33.34;
          input Real Hb;
          Real Hct=Hb/fullHb;

          Real BEp=BE - mHCO3/(1 - Hct) "meq/l";
          Real BEe=BE + mHCO3/Hct "meq/l";

          Real mHCO3;

          Real SID;
          input Real BE;
          input Real pCO2;
          input Real Pi;
          input Real alb;
          output Real pH=fullErythrocyte.pH;
        equation
          plasma.pH = fullErythrocyte.pH;
          BEp = SID - plasma.NSID;
          annotation (experiment(
              StopTime=1,
              __Dymola_NumberOfIntervals=500,
              Tolerance=1e-003));
        end CombinedModelBase;
      end Auxiliary;

      model Wolf_full_blood
        "Implementation of the Wolf acidbase model, reduced to Plasma and Erythrocyte"
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
        parameter Real ZimE=-9.2;
        parameter Real ZimP=-5.3;
        //volumes and CO2 preassure
        parameter Real Vblood(unit="l") = 5;
        input Real pCO2(unit="torr") = 20 + 40*time;
        //
        //derived parameters
        //water volumes in compartments
        parameter Real Vew0(unit="l") = 0.44*0.73*Vblood;
        parameter Real Vpw0(unit="l") = (1 - 0.44)*0.94*Vblood;
        //masses in erythrocytes
        parameter Real m0NaE(unit="mmol") = NaE*Vew0;
        parameter Real m0KE(unit="mmol") = KE*Vew0;
        parameter Real m0ClE(unit="mmol") = ClE*Vew0;
        parameter Real m0Hb(unit="mmol") = Hb*Vew0;
        parameter Real m0DPG(unit="mmol") = DPG*Vew0;
        parameter Real m0ATP(unit="mmol") = ATP*Vew0;
        parameter Real m0GSH(unit="mmol") = GSH*Vew0;
        parameter Real m0imE(unit="mmol") = imE*Vew0;
        parameter Real m0PiE(unit="mmol") = PiE*Vew0;
        //masses in plasma
        parameter Real m0NaP(unit="mmol") = NaP*Vpw0;
        parameter Real m0KP(unit="mmol") = KP*Vpw0;
        parameter Real m0CaP(unit="mmol") = CaP*Vpw0;
        parameter Real m0MgP(unit="mmol") = MgP*Vpw0;
        parameter Real m0ClP(unit="mmol") = ClP*Vpw0;
        parameter Real m0PiP(unit="mmol") = PiP*Vpw0;
        Real m0AlbP(unit="mmol") = AlbP*Vpw0;
        parameter Real m0imP(unit="mmol") = imP*Vpw0;
        //overall masses of mobile ions
        Real MCl(unit="mmol");
        // = m0ClE + m0ClP;
        //
        //masses of mobile ions - 13 unknowns
        Real mClE(unit="mmol");
        Real mClP(unit="mmol",start=m0ClP);
        //
        //volumes of water - 3 unknowns
        Real Vew(unit="l",start=Vew0);
        Real Vpw(unit="l",start=Vpw0);
        //
        // concentrations of bicarbonates - 3 unknowns
        Real HCO3E(unit="mmol/l");
        Real HCO3P(unit="mmol/l");
        //
        //pH dependatn chargers
        Real ZPi=(-1) - 10^(pHP - 6.87)/(1 + 10^(pHP - 6.87));
        Real ZAlb=(-10.7) - 16*(10^(pHP - 7.42)/(1 + 10^(pHP - 7.42)));
        Real ZHb=15.6 - 23*(10^(pHE - 6.69)/(1 + 10^(pHE - 6.69))) - 4*(10^(pHE
             - 7.89)/(1 + 10^(pHE - 7.89))) + 1.5*((1 - 0.75)/0.75);
        Real ZDPG=(-3) - 1*(10^(pHE - 7.56)/(1 + 10^(pHE - 7.56))) - 1*(10^(pHE
             - 7.32)/(1 + 10^(pHE - 7.32)));
        Real ZATP=(-3) - 1*(10^(pHE - 6.8)/(1 + 10^(pHE - 6.8)));
        Real ZGSH=(-1) - 1*(10^(pHE - 8.54)/(1 + 10^(pHE - 8.54))) - 1*(10^(pHE
             - 9.42)/(1 + 10^(pHE - 9.42)));
        Real fiHb=1 + 0.115*C_Hb + 0.0256*C_Hb^2;
        //carbonates
        Real CO3E(unit="mmol/l") = HCO3E*10^(pHE - 10.2);
        Real CO3P(unit="mmol/l") = HCO3P*10^(pHP - 10.2);
        //concentrations
        Real C_NaE(
          unit="mmol/l",
          start=NaE) = m0NaE/Vew;
        Real C_KE(
          unit="mmol/l",
          start=KE) = m0KE/Vew;
        Real C_ClE(
          unit="mmol/l",
          start=ClE) = mClE/Vew;
        Real C_Hb(
          unit="mmol/l",
          start=Hb) = m0Hb/Vew;
        Real C_DPG(
          unit="mmol/l",
          start=DPG) = m0DPG/Vew;
        Real C_ATP(
          unit="mmol/l",
          start=ATP) = m0ATP/Vew;
        Real C_GSH(
          unit="mmol/l",
          start=GSH) = m0GSH/Vew;
        Real C_imE(
          unit="mmol/l",
          start=imE) = m0imE/Vew;
        Real C_PiE(
          unit="mmol/l",
          start=PiE) = m0PiE/Vew;
        //
        Real C_NaP(
          unit="mmol/l",
          start=NaP) = m0NaP/Vpw;
        Real C_KP(
          unit="mmol/l",
          start=KP) = m0KP/Vpw;
        Real C_CaP(
          unit="mmol/l",
          start=CaP) = m0CaP/Vpw;
        Real C_MgP(
          unit="mmol/l",
          start=MgP) = m0MgP/Vpw;
        Real C_ClP(
          unit="mmol/l",
          start=ClP) = mClP/Vpw;
        Real C_PiP(
          unit="mmol/l",
          start=PiP) = m0PiP/Vpw;
        Real C_AlbP(
          unit="mmol/l",
          start=AlbP) = m0AlbP/Vpw;
        Real C_imP(
          unit="mmol/l",
          start=imP) = m0imP/Vpw;
        //osmolality
        Real Oe=(0.93*C_NaE + 0.93*C_KE + 0.93*C_ClE + 0.93*C_PiE + fiHb*C_Hb
             + C_DPG + C_ATP + C_GSH + C_imE + HCO3E + CO3E);
        //*0.73;
        Real Op=(0.93*C_NaP + 0.93*C_KP + 0.93*C_ClP + C_CaP + C_MgP + HCO3P +
            CO3P + 0.93*C_PiP + C_AlbP + C_imP);
        // *0.94;
        //electric charge
        Real Qe=(C_NaE + C_KE - C_ClE - HCO3E - 2*CO3E + ZHb*C_Hb + ZDPG*C_DPG
             + ZATP*C_ATP + ZGSH*C_GSH + ZimE)*Vew;
        Real Qp=(C_NaP + C_KP + 2*C_CaP + 2*C_MgP - C_ClP - HCO3P - 2*CO3P +
            ZPi*C_PiP + ZAlb*C_AlbP + ZimP)*Vpw + X;
        parameter Real X(unit="mEq") = 0;
        //
        input Real BE=0;
        //Real SID;
        Real fH=(Vew + Vew/0.73*(1 - 0.73))/(Vew + Vew/0.73*(1 - 0.73) + Vpw +
            Vpw/0.94*(1 - 0.94));
        Real rCl=C_ClE/C_ClP;
        Real rHCO3=HCO3E/HCO3P;
        //addition of species
        parameter Real XCl(unit="mmol") = 0;

        Real HE(start=10^(-7.2));
        Real HP(start=10^(-7.37));
        //other unknowns
        //pH
        Real pHE(start=7.22) = -log10(HE);
        Real pHP(start=7.4) = -log10(HP);
      equation
        //mass conservation - 6 equations
        MCl = mClE + mClP - XCl;
        //volume conservation - 1 equation
        Vew0 + Vpw0 = Vew + Vpw;
        //
        //donnan equilibrium - 7 equations
        C_ClE/C_ClP = HP/HE;
        //
        //electroneutrality - 3 equations
        Qe = 0;
        Qp = 0;
        //osmotic equilibrium - 2 equations
        Op = Oe;
        //
        //carbonates and pH
        HCO3E = 0.026*pCO2*10^(pHE - 6.11);
        HCO3P = 0.0306*pCO2*10^(pHP - 6.11);
        //
        BE = (1 - 0.023*9)*(HCO3P - 24.4 + (2.3*9 + 7.7)*(pHP - 7.4));
        //SID = (1 - (1 - HCO3E / HCO3P) * fH * fB) * HCO3P + (1 - fH * fB) * (C_AlbP * (8 * pHP - 41) + C_PiP * (0.3 * pHP - 0.4)) + C_Hb * fB * (10.2 * pHP - 73.6) + C_DPG * fH * fB * (0.7 * pHP - 0.5);
        annotation (Icon(coordinateSystem(
              extent={{-100,-100},{100,100}},
              preserveAspectRatio=true,
              initialScale=0.1,
              grid={2,2})), Diagram(coordinateSystem(
              extent={{-100,-100},{100,100}},
              preserveAspectRatio=true,
              initialScale=0.1,
              grid={2,2})));
      end Wolf_full_blood;

      model CompareCombinations
        import FullBloodAcidBase;
        FullBloodAcidBase.FullBloodCombined.comparisson.Auxiliary.CombinedModelBase
          combinedModelBase(
          Hb=15,
          BE=BE,
          pCO2=pCO2,
          Pi=1.15,
          alb=4.4,
          redeclare FullBloodAcidBase.PlasmaElectrochemical.PlasmaFencl plasma)
          annotation (Placement(transformation(extent={{-60,60},{-40,80}})));
        FullBloodAcidBase.FullBloodCombined.comparisson.Auxiliary.CombinedModelBase
          combinedModelBase1(
          Hb=15,
          BE=BE,
          pCO2=pCO2,
          Pi=1.15,
          alb=4.4)
          annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
        FullBloodAcidBase.FullBloodCombined.comparisson.Auxiliary.CombinedModelBase
          combinedModelBase2(
          Hb=15,
          BE=BE,
          pCO2=pCO2,
          Pi=1.15,
          alb=4.4,
          redeclare FullBloodAcidBase.PlasmaElectrochemical.PlasmaWolf plasma)
          annotation (Placement(transformation(extent={{-60,-20},{-40,0}})));

        Real BE=-20 + time*40;
        parameter Real pCO2=40;
      equation

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end CompareCombinations;
    end comparisson;

  end FullBloodCombined;

  package Figures
    extends Modelica.Icons.ExamplesPackage;

    model Figure1_BE_curve "BE curve for graphs with SA nomogram"
      extends Modelica.Icons.Example;
      Real BE=time*60 - 30;
      output Real pCO2BE(start=40) = pco2BECoef[1]*BE^6 + pco2BECoef[2]*BE^5 +
        pco2BECoef[3]*BE^4 + pco2BECoef[4]*BE^3 + pco2BECoef[5]*BE^2 +
        pco2BECoef[6]*BE + pco2BECoef[7];
      output Real pHBE(start=7) = pHBECoef[1]*BE^4 + pHBECoef[2]*BE^3 +
        pHBECoef[3]*BE^2 + pHBECoef[4]*BE + pHBECoef[5];

    protected
      constant Real pco2BECoef[:]={8.3975e-009,-513.9503e-009,3.8105e-006,
          231.6888e-006,-46.5581e-003,353.7105e-003,39.9871e+000};
      constant Real pHBECoef[:]={131.3315e-009,2.5027e-006,175.6144e-006,
          11.9273e-003,7.4001e+000};

      // Dymola plot script

      /*
  createPlot(id=2, position={219, 4, 483, 300}, x="pCO2BE", y={"pHBE"}, range={20.0, 60.00000000000001, 7.1000000000000005, 7.800000000000001}, autoscale=false, grid=true, legend=false, filename="Figure1_4_BE_curve.mat", logX=true, leftTitleType=0, bottomTitleType=0, colors={{28,108,200}}, rightTitleType=0);
  */

    end Figure1_BE_curve;

    model Figure_1
      extends Modelica.Icons.Example;
      Real logpCO2=log10(pCO2);
      Real pCO2=time*60 + 20;
      parameter Real Hb=15;

      FullBloodCombined.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE_10(
        BE=-15,
        pCO2=pCO2,
        Hb=Hb)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
      FullBloodCombined.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE0(
        BE=-0,
        pCO2=pCO2,
        Hb=Hb)
        annotation (Placement(transformation(extent={{-58,20},{-38,40}})));
      FullBloodCombined.comparisson.Auxiliary.ResultSetAtBE resultSetAtBE10(
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

      parameter Real Pi=1.15;
      parameter Real alb=4.4;
      parameter Real Hb=15;
      Real BE=time*40 - 20;
      Real BEfixed=0;
      Real pCO2=time*60 + 20;
      Real pCO2fixed=40;

      FullBloodAcidBase.FullBloodCombined.comparisson.Auxiliary.SetAtAlb
        fixedBE(
        BE=BEfixed,
        pCO2=pCO2,
        alb=alb,
        Hb=Hb,
        Pi=Pi) "Normal model set at fixed BE";
      FullBloodAcidBase.FullBloodCombined.comparisson.Auxiliary.SetAtAlb
        fixedpCO2(
        BE=BE,
        pCO2=pCO2fixed,
        alb=alb,
        Hb=Hb,
        Pi=Pi) "Normal model set at fixed pCO2";

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
      parameter Real pCO2(displayUnit="mmHg") = 40;
      Real BE(displayUnit="meq/l") = ndp.SID*dilutionFactor - normalPlasma.SID
        "Shift in BE for diluted system.";
      Real newSID(displayUnit="meq/l") = ndp.SID*dilutionFactor
        "SID to which the system is diluted to. Only for comparison.";
      // normal SID is 38.97
      parameter Real alb(unit="g/dl") = 4.4;
      parameter Real Pi(unit="mmol/l") = 1.15;
      Real dilutionFactor=1.25 - time;
      Real Hb(unit="g/dl") = 15*dilutionFactor;
      Real HbPerCent=Hb/15*100;
      FullBloodCombined.comparisson.Auxiliary.SetAtAlb dilutedClosed(
        BE=BE,
        pCO2=pCO2*dilutionFactor,
        alb=alb*dilutionFactor,
        Hb=Hb,
        Pi=Pi*dilutionFactor) "Diluted closed system with diluted pCO2"
        annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));
      FullBloodCombined.comparisson.Auxiliary.SetAtAlb dilutedOpenSystem(
        BE=BE,
        pCO2=pCO2,
        alb=alb*dilutionFactor,
        Hb=Hb,
        Pi=Pi*dilutionFactor)
        "Diluted open system, i.e. the pCO2 is held constant by the respiration"
        annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));
      FullBloodCombined.comparisson.Auxiliary.SetAtAlb nondilluted(
        BE=BE,
        pCO2=pCO2,
        alb=alb,
        Hb=Hb,
        Pi=Pi) "Non-dilluted system for comparison"
        annotation (Placement(transformation(extent={{0,-40},{20,-20}})));

    protected
      PlasmaElectrochemical.Auxilliary.FiggeFenclNSID normalPlasma(
        pH0=7.4,
        pCO20=40,
        alb0=alb*dilutionFactor,
        Pi0=Pi*dilutionFactor)
        "normal plasma for establishing NSID in diluted blood";

      PlasmaElectrochemical.Auxilliary.FiggeFenclNSID ndp(
        pH0=7.4,
        pCO20=40,
        alb0=alb,
        Pi0=Pi)
        "normal plasma for establishing NSID in reference, i.e. non-diluted blood";

    public
      PlasmaElectrochemical.Auxilliary.FiggeFencl3 figgeFencl3Open(
        SID=ndp.SID*dilutionFactor,
        pCO2=pCO2,
        Pi=Pi*dilutionFactor,
        alb=alb*dilutionFactor)
        "Comparison to dilution of FiggeFencl3 model in open system, i.e. with constant pCO2"
        annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
      PlasmaElectrochemical.Auxilliary.FiggeFencl3 figgeFencl3Closed(
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

      Real BEf=0 "BE for fixed BE during Alb change";
      Real BEsidf=normalPlasma.SID - albPlasma.SID
        "BE for fixed SID during Alb change";

      parameter Real pCO2(displayUnit="mmHg") = 40;
      parameter Real Hb(unit="g/dl") = 15;
      parameter Real Pi(unit="mmol/l") = 1.15;
      parameter Real alb0(unit="g/dl") = 4.4;
      Real alb(unit="g/dl") = alb0/2 + alb0*time
        "Variable alb from (0.5 - 1.5)*alb0 during 1s";

      FullBloodAcidBase.FullBloodCombined.comparisson.Auxiliary.SetAtAlb
        fixedBE(
        BE=BEf,
        pCO2=pCO2,
        alb=alb,
        Hb=Hb,
        Pi=Pi)
        "Normal model set with fixed BE, that is Albumin is added (or removed) with according H+, resulting in HCO3 increase."
        annotation (Placement(transformation(extent={{20,20},{40,40}})));

      FullBloodAcidBase.FullBloodCombined.comparisson.Auxiliary.SetAtAlb
        fixedSID(
        BE=BEsidf,
        pCO2=pCO2,
        alb=alb,
        Hb=Hb,
        Pi=Pi)
        "Normal model set with fixed SID, that is Albumin is added (or removed) with dissociated cation (e.g. Na, K..)"
        annotation (Placement(transformation(extent={{20,20},{40,40}})));

    protected
      PlasmaElectrochemical.Auxilliary.FiggeFenclNSID albPlasma(
        pH0=7.4,
        pCO20=40,
        alb0=alb,
        Pi0=Pi) "Normal SID for changed albumin level";

      PlasmaElectrochemical.Auxilliary.FiggeFenclNSID normalPlasma(
        pH0=7.4,
        pCO20=40,
        alb0=alb0,
        Pi0=Pi) "Normal standard SID";

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
      Real pCO2=time*60 + 20;
      parameter Real Hb=15;
      parameter Real Hct=Hb/33.34;
      FullBloodEmpirical.Tests.Blood_pCO2 sA_comparison_pCO2_BE_10(
        BEox=-15,
        pCO2=pCO2,
        Hct=Hct)
        annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
      FullBloodEmpirical.Tests.Blood_pCO2 sA_comparison_pCO2_BE0(
        BEox=0,
        pCO2=pCO2,
        Hct=Hct)
        annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
      FullBloodEmpirical.Tests.Blood_pCO2 sA_comparison_pCO2_BE10(
        BEox=15,
        pCO2=pCO2,
        Hct=Hct) annotation (Placement(transformation(extent={{0,20},{20,40}})));
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

      FullBloodEmpirical.SAoriginal plasma(
        Hct=0,
        BEox=BEp,
        pCO2=pCO2)
        annotation (Placement(transformation(extent={{-20,20},{0,40}})));
      FullBloodEmpirical.SAoriginal ery(
        Hct=1,
        BEox=BEe,
        pCO2=pCO2)
        annotation (Placement(transformation(extent={{20,20},{40,40}})));
      FullBloodEmpirical.SAoriginal SABlood(
        Hct=Hct,
        BEox=BE,
        pCO2=pCO2)
        annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
      constant Real fullHb=33.34;
      parameter Real Hct=15/33.34;
      Real BEp(unit="meq/l") = BE - mHCO3/(1 - Hct);
      Real BEe(unit="meq/l") = BE + mHCO3/Hct;
      Real mHCO3;
      //Real pHFullBlood;

      Real BE;
      Real pCO2;

    equation
      plasma.pH = ery.pH;
      BE = 0;
      //time*60 -30;
      //  pCO2 = time*20 + 20;
      pCO2 = time*40 + 20;
      annotation (experiment(
          StopTime=1,
          __Dymola_NumberOfIntervals=500,
          Tolerance=1e-003));
    end testpCO2;

    model Test_Combo_Wolf
      Real pCO2=20 + time*60;
      Real BE=0;
      parameter Real a=0.65/0.94;
      FullBloodCombined.CombinedModel full_blood_combo(pCO2=pCO2, BE=BE)
        annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
      FullBloodCombined.comparisson.Wolf_full_blood wolf_full_blood(
        pCO2=pCO2,
        BE=BE,
        AlbP=a) annotation (Placement(transformation(extent={{0,0},{20,20}})));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Test_Combo_Wolf;

    model Test_Combo_Wolf_at_BE
      Real pCO2=20 + time*60;
      Test_Combo_Wolf test_Combo_Wolf_15(pCO2=pCO2, BE=-15)
        annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
      Test_Combo_Wolf test_Combo_Wolf0(pCO2=pCO2, BE=0)
        annotation (Placement(transformation(extent={{-20,0},{0,20}})));
      Test_Combo_Wolf test_Combo_Wolf15(pCO2=pCO2, BE=15)
        annotation (Placement(transformation(extent={{20,0},{40,20}})));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));

      /*        
// Dymola plot command

createPlot(id=1, position={15, 10, 584, 420}, x="pCO2", y={"test_Combo_Wolf_15.full_blood_combo.pH", "test_Combo_Wolf_15.wolf_full_blood.pHP", "test_Combo_Wolf0.full_blood_combo.pH", "test_Combo_Wolf0.wolf_full_blood.pHP", "test_Combo_Wolf15.wolf_full_blood.pHP", "test_Combo_Wolf15.full_blood_combo.pH"}, range={20.0, 80.0, 7.0, 7.8}, autoscale=false, grid=true, legend=false, filename="dsres.mat", logX=true, leftTitleType=2, leftTitle="pH", colors={{238,46,47}, {0,0,0}, {238,46,47}, {0,0,0}, {0,0,0}, {238,46,47}}, patterns={LinePattern.Solid, LinePattern.Dash, LinePattern.Solid, LinePattern.Dash, LinePattern.Dash, LinePattern.Solid}, thicknesses={0.5, 0.5, 0.5, 0.5, 0.5, 0.5});
*/
    end Test_Combo_Wolf_at_BE;
  end Tests;

  package FullBloodSubmodelComponent
    extends Modelica.Icons.ExamplesPackage;
    package Components
      model CombinedModel
        constant Real mmHgToPa(unit="1/Pa") = 1/133;
        constant Real albConversion=1000*10/66500;
        constant Real HbConversion=16000/1000/10
          "Conversion from mol/m3 of single hem to g/dl. One hem has around 16 kDa.";
        FullBloodCombined.CombinedModel combinedModel(
          BE=BEox,
          pCO2=pCO2*mmHgToPa,
          alb=Alb*albConversion,
          Hb=ctHb*HbConversion,
          Pi=Pi) annotation (Placement(transformation(extent={{-20,0},{0,20}})));
        Physiolibrary.Types.RealIO.PressureInput pCO2
          annotation (Placement(transformation(extent={{-120,30},{-80,70}})));
        Physiolibrary.Types.RealIO.ConcentrationInput BEox
          annotation (Placement(transformation(extent={{-120,-10},{-80,30}})));
        Physiolibrary.Types.RealIO.ConcentrationInput Alb
          annotation (Placement(transformation(extent={{-120,-50},{-80,-10}})));
        Physiolibrary.Types.RealIO.ConcentrationInput Pi
          annotation (Placement(transformation(extent={{-120,-90},{-80,-50}})));
        Physiolibrary.Types.RealIO.pHOutput pH=combinedModel.pH;
        Physiolibrary.Types.RealIO.ConcentrationInput ctHb
          annotation (Placement(transformation(extent={{-40,70},{0,110}})));
        annotation (
          Placement(transformation(extent={{80,0},{100,20}})),
          Documentation(info="<html>
<p><span style=\"font-family: Arial,sans-serif; color: #222222; background-color: #ffffff;\">Interface using the SI units, as the Combined model relies on units originally  used by the cited models.</span></p>
</html>"),
          Icon(coordinateSystem(preserveAspectRatio=false)),
          Diagram(coordinateSystem(preserveAspectRatio=false)));
      end CombinedModel;

      package OSA

        function aFrom
          input Real pH;
          input Real pCO2;
          input Real MetHb;
          input Real HbF;
          input Real cDPG;
          output Real returnValue;
        protected
          Real dadpH=-0.88;
          Real dadlnpCO2=0.048;
          Real dadxMetHb=-0.7;
          Real dadxHbF=-0.25;
          Real dadcDPG0=0.3;
          Real pH0=7.40;
          Real pCO20=5.33;
          Real dadcDPGxHbF=-0.1;
          Real cDPG0=5.0;
        algorithm
          returnValue := dadpH*(pH - pH0) + dadlnpCO2*log(pCO2/pCO20) +
            dadxMetHb*MetHb + dadxHbF*HbF + (dadcDPG0 + dadcDPGxHbF*HbF)*(cDPG/
            cDPG0 - 1.0);
        end aFrom;

        function sCO
          input Real FCOHb;
          input Real FMetHb;
          output Real returnValue;
        protected
          Real xFCOHb;
        algorithm
          if FCOHb < 0 then
            xFCOHb := 0;
          else
            xFCOHb := FCOHb;
          end if;
          returnValue := xFCOHb/(1.0 - FMetHb);
        end sCO;

        function antilogit
          input Real x;
          output Real returnValue;
        algorithm
          returnValue := exp(x)/(1.0 + exp(x));
        end antilogit;

        function logit
          input Real x;
          output Real returnValue;
        algorithm
          returnValue := log(x/(1 - x));
        end logit;

        function h
          input Real a;
          output Real returnValue;
        protected
          Real h0=3.5;
        algorithm
          returnValue := h0 + a;
        end h;

        function x
          input Real pO2CO;
          input Real a;
          input Real T;
          output Real returnValue;
        protected
          Real p0=7.0;
          Real T0=37.0;
          Real dbdT=0.055;
        algorithm
          returnValue := log(pO2CO/p0) - a - dbdT*(T - T0);
        end x;

        function dydx
          input Real pO2CO;
          input Real a;
          input Real T;
          output Real returnValue;
        protected
          Real k=0.5342857;
        algorithm
          returnValue := 1 + h(a)*k*(1 - (tanh(k*x(
                    pO2CO,
                    a,
                    T)))^2);
        end dydx;

        function y
          input Real pO2CO;
          input Real a;
          input Real T;
          output Real returnValue;
        protected
          Real y0=1.8747;
          Real k=0.5342857;
        algorithm
          returnValue := y0 + x(
                    pO2CO,
                    a,
                    T) + h(a)*tanh(k*x(
                    pO2CO,
                    a,
                    T));
        end y;

        function sO2CO
          input Real pO2CO;
          input Real a;
          input Real T;
          output Real returnValue;
        algorithm
          returnValue := antilogit(y(
                    pO2CO,
                    a,
                    T));
        end sO2CO;

        function MpCOof
          input Real pO2CO;
          input Real a;
          input Real T;
          input Real FCOHb;
          input Real FMetHb;
          output Real returnValue;
        algorithm
          returnValue := (pO2CO/sO2CO(
                    pO2CO,
                    a,
                    T))*sCO(FCOHb, FMetHb);
        end MpCOof;

        function pO2fr
          input Real sO2;
          input Real a;
          input Real T;
          input Real FCOHb;
          input Real FMetHb;
          output Real returnValue;
        protected
          Real pO2CO;
          Real sO2CO;
          Real ym;
          Real yc;
          Real dydxc;
          Real p0=7.0;
          Real dbdT=0.055;
          Real T0=37;
          Boolean doit;
          Real Epsilon=0.000001;
        algorithm
          pO2CO := exp(log(p0) + a + dbdT*(T - T0));
          sO2CO := sO2 + sCO(FCOHb, FMetHb)*(1 - sO2);
          ym := logit(sO2CO);
          doit := false;
          while not doit loop
            yc := y(  pO2CO,
                      a,
                      T);
            dydxc := dydx(
                      pO2CO,
                      a,
                      T);
            pO2CO := exp(log(pO2CO) + (ym - yc)/dydxc);
            doit := abs(ym - yc) < Epsilon;
          end while;
          returnValue := pO2CO - MpCOof(
                    pO2CO,
                    a,
                    T,
                    FCOHb,
                    FMetHb);
        end pO2fr;

        function sO2fr
          input Real pO2CO;
          input Real a;
          input Real T;
          input Real FCOHb;
          input Real FMetHb;
          output Real returnValue;
        protected
          Real sO2COc;
          Real sCOc;
        algorithm
          sO2COc := sO2CO(
                    pO2CO,
                    a,
                    T);
          sCOc := sCO(FCOHb, FMetHb);
          returnValue := (sO2COc - sCOc)/(1 - sCOc);
        end sO2fr;

        function sO2of "calculation of oxygen hemoglobin saturation"
          input Real pO2T "Po2 at given temperature in kPa";
          input Real pHT "pH at given temperature";
          input Real pCO2T "pCO2 at given temperature in kPa";
          input Real cDPG "2'3 DPG koncentration in mmol/l";
          input Real FCOHb "substance fraction of carboxyhemoglobin";
          input Real FMetHb "substance fraction of hemiglobin";
          input Real FHbF "substance fraction of fetal hemoglobin";
          input Real TPt "temperature in?C";
          output Real returnValue "oxygen hemoglobin saturation";
        protected
          Real MpCOa;
          Real MpCOb;
          Real sCOc;
          Boolean doit;
          Real a;
          Real Epsilon=0.000001;
        algorithm
          a := aFrom(
                    pHT,
                    pCO2T,
                    FMetHb,
                    FHbF,
                    cDPG);
          sCOc := sCO(FCOHb, FMetHb);
          if sCOc > 0 then
            MpCOa := pO2fr(
                      sCOc,
                      a,
                      TPt,
                      0,
                      FMetHb);
          else
            MpCOa := 0;
          end if;
          MpCOb := MpCOa;
          doit := false;
          while (not doit) loop
            MpCOb := 0.6*MpCOa + 0.4*MpCOb;
            MpCOa := MpCOof(
                      pO2T + MpCOb,
                      a,
                      TPt,
                      FCOHb,
                      FMetHb);
            doit := (abs(MpCOa - MpCOb) < Epsilon);
          end while;
          returnValue := sO2fr(
                    pO2T + MpCOa,
                    a,
                    TPt,
                    FCOHb,
                    FMetHb);
        end sO2of;

        function aO2
          input Real temp;
          output Real returnValue;
        algorithm
          returnValue := exp(log(0.0105) - 0.0115*(temp - 37.0) + 0.5*0.00042*(
            temp - 37.0)^2);
        end aO2;

        function dissO2 "concentration of dissolved oxygen in blood"
          input Real pO2;
          input Real temp;
          output Real returnValue "dissolved blood oxygen in mmol/l";
        algorithm
          returnValue := aO2(temp)*pO2;
        end dissO2;

        function ceHbof "effective hemoglobin concentration in mmol/l"
          input Real ctHb "concentration of hemoglobin in mmol/l";
          input Real FCOHb "substance fraction of carboxyhemoglobin";
          input Real FMetHb "substance fraction of hemoglobin";
          output Real returnValue "effective contentration of hemoglobin";
        algorithm
          returnValue := ctHb*(1 - FCOHb - FMetHb);
        end ceHbof;

        function O2totalSI "Calculation of concentration of total oxygen"
          input Real ctHb "conentration of hemoglobin in mmol/l";
          input Real pO2 "pO2 at givent temperature in Pa";
          input Real pHp "pH in plasma at given temperature";
          input Real pCO2 "pCO2 at given temperature in Pa";
          input Real cDPG "concentration of 2,3 diphosphoglycerate in mmol/l";
          input Real FCOHb "substance fraction of carboxyhemoglobin";
          input Real FMetHb "substance fraction of hemiglobin";
          input Real FHbF "substance fraction of fetal hemogobin";
          input Real temp "temperature in °K";
          output Real ctO2
            "concentration of total blood oxygen concentration in mmol/l";
          output Real sO2t
            "oxygen saturation of hemoglobin at given temperature";
          output Real dissO2t
            "koncentration of dissolved oxygen in blood in mmol/l";
          output Real ceHb "effective hemoglobin concentration in mmol/l";
        algorithm

          sO2t := sO2of(
                    pO2/1000,
                    pHp,
                    pCO2/1000,
                    cDPG,
                    FCOHb,
                    FMetHb,
                    FHbF,
                    temp - 273.15);
          ceHb := ceHbof(
                    ctHb,
                    FCOHb,
                    FMetHb);
          dissO2t := dissO2(pO2/1000, temp - 273.15);
          ctO2 := dissO2t + sO2t*ceHb;
        end O2totalSI;

        function lg
          input Real x;
          output Real result;
        algorithm
          result := (log(x))/log(10);
          //it is not necessary, in Modelica exists embeded function log10
        end lg;

        function antilg
          input Real x;
          output Real result;
        algorithm
          result := exp(log(10)*x);
        end antilg;

        function aCO2of
          input Real T;
          output Real result;
        protected
          Real aCO2T0=0.23;
          //mM/kPa
          Real dlgaCO2dT=-0.0092;
          // lg(mM/kPa)/K
          Real T0=37;
        algorithm
          result := aCO2T0*antilg(dlgaCO2dT*(T - T0));
        end aCO2of;

        function pKof
          input Real T;
          output Real result;
        protected
          Real pKT0=6.1;
          Real dpKdT=-0.0026;
          Real T0=37;
        algorithm
          result := pKT0 + dpKdT*(T - T0);
        end pKof;

        function cHCO3of "calculation of plasma bicarbonate concentration"
          input Real pH "plasma pH at given temperature in mmol/l";
          input Real pCO2 "pCO2 in kPa";
          input Real T "temperature in °C";
          output Real HCO3p "plasma bicarbonate concentration in mmol/l";
        algorithm
          HCO3p := pCO2*aCO2of(T)*antilg(pH - pKof(T));
        end cHCO3of;

        function pCO22of
          input Real pCO21;
          input Real T1;
          input Real T2;
          input Real cHb;
          input Real cAlb;
          input Real pH1;
          output Real result;
        protected
          Real betaX;
          Real dpHdT1;
          Real pH2;
          Real cHCO3;
          Real dlgpCO2dT1;
          Real pCO22;
          Real dpHdT2;
          Real dlgpCO2dT2;
          Real dpHdTmean;
          Real dlgpCO2dTmean;
          Real cAlbN=0.66;
        algorithm
          betaX := 7.7 + 8*(cAlb - cAlbN) + 2.3*cHb;

          dpHdT1 := (-0.0026 - betaX*0.016*(1/(2.3*cHCO3of(
                    pH1,
                    pCO21,
                    T1)) + 1/(2.3*pCO21*aCO2of(T1))))/(1 + betaX*(1/(2.3*
            cHCO3of(pH1,
                    pCO21,
                    T1)) + 1/(2.3*pCO21*aCO2of(T1))));
          pH2 := pH1 + dpHdT1*(T2 - T1);

          cHCO3 := cHCO3of(
                    pH1,
                    pCO21,
                    T1);
          dlgpCO2dT1 := -0.0026 - (-0.0092) - dpHdT1 + (1/(2.3*cHCO3))*(betaX*
            dpHdT1 + betaX*0.016);
          pCO22 := antilg(lg(pCO21) + dlgpCO2dT1*(T2 - T1));

          dpHdT2 := (-0.0026 - betaX*0.016*(1/(2.3*cHCO3of(
                    pH2,
                    pCO22,
                    T2)) + 1/(2.3*pCO22*aCO2of(T2))))/(1 + betaX*(1/(2.3*
            cHCO3of(pH2,
                    pCO22,
                    T2)) + 1/(2.3*pCO22*aCO2of(T2))));
          dpHdTmean := (dpHdT1 + dpHdT2)/2;
          pH2 := pH1 + dpHdTmean*(T2 - T1);

          cHCO3 := cHCO3of(
                    pH2,
                    pCO22,
                    T2);
          dlgpCO2dT2 := -0.0026 - (-0.0092) - dpHdT2 + (1/(2.3*cHCO3))*(betaX*
            dpHdT2 + betaX*0.016);
          dlgpCO2dTmean := (dlgpCO2dT1 + dlgpCO2dT2)/2;

          result := antilg(lg(pCO21) + dlgpCO2dTmean*(T2 - T1));
        end pCO22of;

        function pH2of
          input Real pH1;
          input Real T1;
          input Real T2;
          input Real cHb;
          input Real cAlb;
          input Real pCO21;
          output Real result;
        protected
          Real betaX;
          Real dpHdT1;
          Real pH2;
          Real cHCO3;
          Real dlgpCO2dT1;
          Real pCO22;
          Real dpHdT2;
          Real dpHdTmean;
          Real cAlbN=0.66;
        algorithm
          betaX := 7.7 + 8*(cAlb - cAlbN) + 2.3*cHb;
          dpHdT1 := (-0.0026 - betaX*0.016*(1/(2.3*cHCO3of(
                    pH1,
                    pCO21,
                    T1)) + 1/(2.3*pCO21*aCO2of(T1))))/(1 + betaX*(1/(2.3*
            cHCO3of(pH1,
                    pCO21,
                    T1)) + 1/(2.3*pCO21*aCO2of(T1))));
          pH2 := pH1 + dpHdT1*(T2 - T1);

          cHCO3 := cHCO3of(
                    pH1,
                    pCO21,
                    T1);
          dlgpCO2dT1 := -0.0026 - (-0.0092) - dpHdT1 + (1/(2.3*cHCO3))*(betaX*
            dpHdT1 + betaX*0.016);
          pCO22 := antilg(lg(pCO21) + dlgpCO2dT1*(T2 - T1));

          dpHdT2 := (-0.0026 - betaX*0.016*(1/(2.3*cHCO3of(
                    pH2,
                    pCO22,
                    T2)) + 1/(2.3*pCO22*aCO2of(T2))))/(1 + betaX*(1/(2.3*
            cHCO3of(pH2,
                    pCO22,
                    T2)) + 1/(2.3*pCO22*aCO2of(T2))));
          dpHdTmean := (dpHdT1 + dpHdT2)/2;
          result := pH1 + dpHdTmean*(T2 - T1);

        end pH2of;

        function CO2totalSI "Calculation of blood total CO2 concentration"
          input Real pH "plasma pH at given temperature";
          input Real pCO2 "pCO2 at given temperatura in Pa";
          input Real T "temperature in ?C";
          input Real ctHb "Hemoglobin concentration in mmol/l";
          input Real sO2 "O2 hemoglobin saturation";
          output Real ctCO2B "Total blood CO2 concetratoin in mmol/l";
          output Real cHCO3 "plasma concentration of bicarbonate in mmol/l";
          output Real dCO2 "dissolved CO2 concentration in plasma";
        protected
          Real dpHEdpHP=0.77;
          Real dpHEdsO2=0.035;
          Real pHEx=7.84;
          Real sO2x=0.06;
          Real aCO2E0=0.195;
          Real ctHbE=21;
          Real pHE0=7.19;
          Real pKE0=6.125;
          Real pHT0;
          Real pCO2T0;
          Real pKE;
          Real pHE;
          Real ctCO2E;
          Real phiEB;
          Real T0=37;
          Real cAlbN=0.66;
          Real cAlb;
          Real pH0=7.40;
          Real aCO2;
          Real tCO2p;

        algorithm
          // pCO2T0 := pCO22of (pCO2, T, T0, ctHb);
          cAlb := cAlbN;
          // albumin has minimal influence on total CO2 concentration
          pCO2T0 := pCO22of(
                    pCO2/1000,
                    T - 273.15,
                    T0,
                    ctHb,
                    cAlb,
                    pH);
          // pHT0 := pH2of (pH, T, T0, ctHb);
          pHT0 := pH2of(
                    pH,
                    T - 273.15,
                    T0,
                    ctHb,
                    cAlb,
                    pCO2);
          pHE := pHE0 + dpHEdpHP*(pHT0 - pH0) + dpHEdsO2*(1 - sO2);
          //or : (pHE - 6.9) = alpha*(pHP - pH0), where alpha = 0.7 + f*(1 - sO2)
          pKE := pKE0 - lg(1 + antilg(pHE - pHEx - sO2x*sO2));
          ctCO2E := aCO2E0*pCO2T0*(1 + antilg(pHE - pKE));
          phiEB := ctHb/ctHbE;
          // !! !! it is hematokrit!!!!!!!
          //tCO2p := pCO2T0 * aCO2of(T0)*(1 + antilg(pHT0-pKof(T0)));
          aCO2 := aCO2of(T0);
          cHCO3 := aCO2*pCO2T0*antilg(pHT0 - pKof(T0));
          dCO2 := aCO2*pCO2T0;
          ctCO2B := ctCO2E*phiEB + (dCO2 + cHCO3)*(1 - phiEB);
          //ctCO2B :=  ctCO2E*phiEB + ctCO2Pof( pHT0, pCO2T0, T0)*(1 - phiEB);
          //ctCO2B :=  ctCO2E*phiEB + tCO2p*(1 - phiEB);
        end CO2totalSI;

        model ctO2content

          Physiolibrary.Types.RealIO.pHInput pH annotation (Placement(
                transformation(extent={{-120,70},{-80,110}}),
                iconTransformation(extent={{-120,32},{-100,52}})));
          Physiolibrary.Types.RealIO.PressureInput pCO2(start=5330) annotation
            (Placement(transformation(extent={{-120,20},{-80,60}}),
                iconTransformation(extent={{-120,-10},{-100,10}})));
          Physiolibrary.Types.RealIO.TemperatureInput T(start=310.15)
            annotation (Placement(transformation(extent={{-120,-20},{-80,20}}),
                iconTransformation(extent={{-120,-50},{-100,-30}})));
          Physiolibrary.Types.RealIO.PressureInput pO2 annotation (Placement(
                transformation(extent={{-132,-54},{-92,-14}}),
                iconTransformation(extent={{-120,70},{-100,90}})));
          Physiolibrary.Types.RealIO.FractionInput FCOHb annotation (Placement(
                transformation(extent={{60,-100},{100,-60}}),
                iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=180,
                origin={110,-80})));
          Physiolibrary.Types.RealIO.FractionInput FHbF annotation (Placement(
                transformation(extent={{60,-60},{100,-20}}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=180,
                origin={110,-40})));
          Physiolibrary.Types.RealIO.FractionInput FMetHb annotation (Placement(
                transformation(extent={{60,-20},{100,20}}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=180,
                origin={110,0})));
          Physiolibrary.Types.RealIO.ConcentrationInput cDPG annotation (
              Placement(transformation(extent={{60,20},{100,60}}),
                iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=180,
                origin={110,40})));
          Physiolibrary.Types.RealIO.ConcentrationInput ctHb annotation (
              Placement(transformation(extent={{60,60},{100,100}}),
                iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=180,
                origin={110,80})));

          Physiolibrary.Types.RealIO.FractionOutput sO2 annotation (Placement(
                transformation(extent={{-30,-112},{10,-72}}),
                iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={60,-110})));
          Physiolibrary.Types.RealIO.ConcentrationOutput totalO2 annotation (
              Placement(transformation(extent={{-80,-100},{-60,-80}}),
                iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={20,-110})));
          Physiolibrary.Types.RealIO.ConcentrationOutput cdO2p
            "dissolved O2 concentration in plasma" annotation (Placement(
                transformation(extent={{-80,-100},{-60,-80}}),
                iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-20,-110})));
          Physiolibrary.Types.RealIO.ConcentrationOutput ceHb
            "effective concentration of hemoglobin" annotation (Placement(
                transformation(extent={{-80,-100},{-60,-80}}),
                iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-60,-110})));

        algorithm
          (totalO2,sO2,cdO2p,ceHb) := O2totalSI(
                    ctHb,
                    pO2,
                    pH,
                    pCO2,
                    cDPG,
                    FCOHb,
                    FMetHb,
                    FHbF,
                    T);

          annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={
                    {-100,-100},{100,100}}),graphics={Rectangle(
                          extent={{-100,100},{100,-100}},
                          lineColor={28,108,200},
                          fillColor={255,255,0},
                          fillPattern=FillPattern.Solid),Text(
                          extent={{-60,66},{64,-34}},
                          lineColor={28,108,200},
                          textString="O2 total")}), Diagram(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}})));
        end ctO2content;

        model ctCO2content

          Physiolibrary.Types.RealIO.PressureInput pCO2(start=5330)
            "pCO2 in Pa" annotation (Placement(transformation(extent={{-120,20},
                    {-80,60}}), iconTransformation(extent={{-120,70},{-100,90}})));
          Physiolibrary.Types.RealIO.pHInput pH annotation (Placement(
                transformation(extent={{-120,70},{-80,110}}),
                iconTransformation(extent={{-120,30},{-100,50}})));
          Physiolibrary.Types.RealIO.TemperatureInput T(start=310.15)
            "temperature (in Kelvins)" annotation (Placement(transformation(
                  extent={{-120,-20},{-80,20}}), iconTransformation(extent={{-120,
                    -10},{-100,10}})));
          Physiolibrary.Types.RealIO.ConcentrationInput ctHb
            "hemoglobin concentration (mmol/l)" annotation (Placement(
                transformation(extent={{60,60},{100,100}}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-110,-40})));
          Physiolibrary.Types.RealIO.FractionInput sO2
            "O2 hemoglobin saturation " annotation (Placement(transformation(
                  extent={{-120,-98},{-80,-58}}), iconTransformation(extent={{-120,
                    -90},{-100,-70}})));
          Physiolibrary.Types.RealIO.ConcentrationOutput ctCO2
            "total blood CO2 concentration (in mmol/l)" annotation (Placement(
                transformation(extent={{100,30},{120,50}}), iconTransformation(
                  extent={{100,30},{120,50}})));
          Physiolibrary.Types.RealIO.ConcentrationOutput cHCO3
            "plasma HCO3 concentration (in mmol/l)" annotation (Placement(
                transformation(extent={{100,60},{120,80}}), iconTransformation(
                  extent={{100,-10},{120,10}})));
          Physiolibrary.Types.RealIO.ConcentrationOutput cdCO2p
            "plasma CO2 dissolved concentration (in mmol/l)" annotation (
              Placement(transformation(extent={{100,60},{120,80}}),
                iconTransformation(extent={{100,-50},{120,-30}})));

        algorithm
          (ctCO2,cHCO3,cdCO2p) := CO2totalSI(
                    pH,
                    pCO2,
                    T,
                    ctHb,
                    sO2);

          annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={
                    {-100,-100},{100,100}}),graphics={Rectangle(
                          extent={{-100,100},{100,-100}},
                          lineColor={28,108,200},
                          fillColor={255,255,0},
                          fillPattern=FillPattern.Solid),Text(
                          extent={{-66,46},{82,-26}},
                          lineColor={28,108,200},
                          fillColor={255,255,0},
                          fillPattern=FillPattern.Solid,
                          textString="CO2 total")}));
        end ctCO2content;
        annotation (Documentation(info="<html>
<p>Modelica algorithm reiplmementation of Siggaard-Andersen Oxygen Status Algorithm, based on original published source codes.</p>
<p><br>Sources:</p>
<p><span style=\"font-family: Arial,sans-serif; color: #222222; background-color: #ffffff;\">SIGGAARD‐ANDERSEN, M. A. D. S., and Ole SIGGAARD‐ANDERSEN. &quot;Oxygen status algorithm, version 3, with some applications.&quot;&nbsp;<i>Acta Anaesthesiologica Scandinavica</i>&nbsp;39.s107 (1995): 13-20.</span></p>
</html>"));
      end OSA;
    end Components;

    model FullBlood
      Physiolibrary.Chemical.Components.Substance totalO2(useNormalizedVolume=
            false, solute_start=0.008253)
        annotation (Placement(transformation(extent={{-82,60},{-62,80}})));
      Physiolibrary.Chemical.Components.Substance totalCO2(useNormalizedVolume=
            false, solute_start=0.021646)
        annotation (Placement(transformation(extent={{-82,20},{-62,40}})));
      Physiolibrary.Chemical.Components.Substance totalBE(useNormalizedVolume=
            false, solute_start=0)
        annotation (Placement(transformation(extent={{-82,-20},{-62,0}})));
      Physiolibrary.Chemical.Interfaces.ChemicalPort_a port_O2
        annotation (Placement(transformation(extent={{-100,60},{-80,80}})));
      Physiolibrary.Chemical.Interfaces.ChemicalPort_a port_CO2
        annotation (Placement(transformation(extent={{-100,20},{-80,40}})));
      Physiolibrary.Chemical.Interfaces.ChemicalPort_a port_BE
        annotation (Placement(transformation(extent={{-100,-20},{-80,0}})));
      Components.OSA.ctO2content ctO2content
        annotation (Placement(transformation(extent={{40,40},{80,80}})));
      Components.OSA.ctCO2content ctCO2content
        annotation (Placement(transformation(extent={{40,-20},{80,20}})));
      Components.CombinedModel combinedModel
        annotation (Placement(transformation(extent={{40,-80},{80,-40}})));

      Physiolibrary.Chemical.Sensors.ConcentrationMeasure concentrationMeasure
        annotation (Placement(transformation(extent={{-60,60},{-40,80}})));
      Physiolibrary.Chemical.Sensors.ConcentrationMeasure concentrationMeasure1
        annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
      Physiolibrary.Chemical.Sensors.ConcentrationMeasure concentrationMeasure2
        annotation (Placement(transformation(extent={{-60,-20},{-40,0}})));
      Modelica.Blocks.Math.InverseBlockConstraints inverseBlockConstraints
        annotation (Placement(transformation(extent={{-26,54},{-2,70}})));
      Modelica.Blocks.Math.InverseBlockConstraints inverseBlockConstraints1
        annotation (Placement(transformation(extent={{-26,14},{-2,30}})));
      Physiolibrary.Types.Constants.ConcentrationConst cDPG(k=5)
        annotation (Placement(transformation(extent={{96,64},{87,70}})));
      Physiolibrary.Types.Constants.FractionConst fMetHb(k=0.005)
        annotation (Placement(transformation(extent={{96,56},{87,62}})));
      Physiolibrary.Types.Constants.FractionConst fHbF(k=0.005)
        annotation (Placement(transformation(extent={{96,48},{87,54}})));
      Physiolibrary.Types.Constants.FractionConst fCOHb(k=0.005)
        annotation (Placement(transformation(extent={{96,40},{87,46}})));


      Physiolibrary.Types.RealIO.VolumeInput volume if UseVolInput
        "total compartment volume" annotation (Placement(transformation(extent=
                {{60,60},{100,100}}),iconTransformation(extent={{60,60},{100,
                100}})));
      Physiolibrary.Types.RealIO.TemperatureInput T if UseTInput annotation (
          Placement(transformation(extent={{-16,-14},{12,14}}),
            iconTransformation(extent={{60,20},{100,60}})));
      Physiolibrary.Types.RealIO.ConcentrationInput ctHb if UsectHctInput
        "Total concentraion of hemoglobin hem (i.e. 4 times the hemoglobin)"
        annotation (Placement(transformation(extent={{-16,-40},{14,-10}}),
            iconTransformation(extent={{60,-20},{100,20}})));
      Physiolibrary.Types.RealIO.ConcentrationInput Alb if UseAlbInput
        annotation (Placement(transformation(extent={{-16,-82},{14,-52}}),
            iconTransformation(extent={{60,-60},{100,-20}})));
      Physiolibrary.Types.RealIO.ConcentrationInput Pi if UsePiInput
        annotation (Placement(transformation(extent={{-16,-98},{14,-68}}),
            iconTransformation(extent={{60,-100},{100,-60}})));

      parameter Physiolibrary.Types.Volume Vol0=0.001;
      parameter Physiolibrary.Types.Temperature T0=310.15;
      parameter Physiolibrary.Types.Concentration ctHb0=60;
      parameter Physiolibrary.Types.Concentration Alb0=0.66;
      parameter Physiolibrary.Types.Concentration Pi0=1.11;
      parameter Boolean UseVolInput=true;
      parameter Boolean UseTInput=true;
      parameter Boolean UsectHctInput=true;
      parameter Boolean UseAlbInput=true;
      parameter Boolean UsePiInput=true;
    equation
      if not UseVolInput then
        totalO2.solutionVolume = Vol0;
        totalCO2.solutionVolume = Vol0;
        totalBE.solutionVolume = Vol0;
      end if;
      if not UseTInput then
        ctO2content.T = T0;
        ctCO2content.T = T0;
      end if;
      if not UsectHctInput then
        combinedModel.ctHb = ctHb0;
        ctO2content.ctHb = ctHb0;
        ctCO2content.ctHb = ctHb0;
      end if;
      if not UseAlbInput then
        combinedModel.Alb = Alb0;
      end if;
      if not UsePiInput then
        combinedModel.Pi = Pi0;
      end if;

      connect(port_O2, totalO2.q_out) annotation (Line(
          points={{-90,70},{-72,70}},
          color={107,45,134},
          thickness=1));
      connect(port_CO2, totalCO2.q_out) annotation (Line(
          points={{-90,30},{-72,30}},
          color={107,45,134},
          thickness=1));
      connect(port_BE, totalBE.q_out) annotation (Line(
          points={{-90,-10},{-72,-10}},
          color={107,45,134},
          thickness=1));
      connect(totalO2.solutionVolume, volume)
        annotation (Line(points={{-76,74},{-76,80},{80,80}}, color={0,0,127}));
      connect(totalCO2.solutionVolume, volume)
        annotation (Line(points={{-76,34},{-76,80},{80,80}}, color={0,0,127}));
      connect(totalBE.solutionVolume, volume)
        annotation (Line(points={{-76,-6},{-76,80},{80,80}}, color={0,0,127}));
      connect(combinedModel.pH, ctO2content.pH) annotation (Line(points={{78,-58},
              {100,-58},{100,96},{20,96},{20,68.4},{38,68.4}}, color={0,0,127}));
      connect(combinedModel.pH, ctCO2content.pH) annotation (Line(points={{78,-58},
              {100,-58},{100,96},{20,96},{20,8},{38,8}},color={0,0,127}));
      connect(ctO2content.sO2, ctCO2content.sO2) annotation (Line(points={{72,
              38},{72,28},{24,28},{24,-16},{38,-16}}, color={0,0,127}));
      connect(T, ctCO2content.T) annotation (Line(points={{-2,1.77636e-015},{18,
              1.77636e-015},{18,0},{38,0}}, color={0,0,127}));
      connect(T, ctO2content.T) annotation (Line(points={{-2,1.77636e-015},{14,
              1.77636e-015},{14,52},{38,52}}, color={0,0,127}));
      connect(ctHb, combinedModel.ctHb) annotation (Line(points={{-1,-25},{56,-25},
              {56,-42}}, color={0,0,127}));
      connect(ctHb, ctCO2content.ctHb) annotation (Line(points={{-1,-25},{32,-25},
              {32,-8},{38,-8}}, color={0,0,127}));
      connect(ctO2content.ctHb, ctHb) annotation (Line(points={{82,76},{98,76},
              {98,-25},{-1,-25}},color={0,0,127}));
      connect(totalO2.q_out, concentrationMeasure.q_in) annotation (Line(
          points={{-72,70},{-50,70}},
          color={107,45,134},
          thickness=1));
      connect(totalCO2.q_out, concentrationMeasure1.q_in) annotation (Line(
          points={{-72,30},{-50,30}},
          color={107,45,134},
          thickness=1));
      connect(totalBE.q_out, concentrationMeasure2.q_in) annotation (Line(
          points={{-72,-10},{-50,-10}},
          color={107,45,134},
          thickness=1));
      connect(concentrationMeasure2.concentration, combinedModel.BEox)
        annotation (Line(points={{-50,-18},{-50,-58},{40,-58}}, color={0,0,127}));
      connect(concentrationMeasure.concentration, inverseBlockConstraints.u1)
        annotation (Line(points={{-50,62},{-27.2,62}}, color={0,0,127}));
      connect(ctO2content.totalO2, inverseBlockConstraints.u2) annotation (Line(
            points={{64,38},{64,32},{-23.6,32},{-23.6,62}}, color={0,0,127}));
      connect(concentrationMeasure1.concentration, inverseBlockConstraints1.u1)
        annotation (Line(points={{-50,22},{-27.2,22}}, color={0,0,127}));
      connect(ctCO2content.ctCO2, inverseBlockConstraints1.u2) annotation (Line(
            points={{82,8},{84,8},{84,24},{-23.6,24},{-23.6,22}}, color={0,0,
              127}));
      connect(ctO2content.pCO2, inverseBlockConstraints1.y2) annotation (Line(
            points={{38,60},{30,60},{30,22},{-3.8,22}}, color={0,0,127}));
      connect(ctCO2content.pCO2, inverseBlockConstraints1.y2) annotation (Line(
            points={{38,16},{30,16},{30,22},{-3.8,22}}, color={0,0,127}));
      connect(combinedModel.pCO2, inverseBlockConstraints1.y2) annotation (Line(
            points={{40,-50},{30,-50},{30,22},{-3.8,22}}, color={0,0,127}));
      connect(ctO2content.cDPG, cDPG.y) annotation (Line(points={{82,68},{84,68},
              {84,67},{85.875,67}}, color={0,0,127}));
      connect(ctO2content.FMetHb, fMetHb.y) annotation (Line(points={{82,60},{
              84,60},{84,59},{85.875,59}}, color={0,0,127}));
      connect(ctO2content.FHbF, fHbF.y) annotation (Line(points={{82,52},{84,52},
              {84,51},{85.875,51}}, color={0,0,127}));
      connect(ctO2content.FCOHb, fCOHb.y) annotation (Line(points={{82,44},{84,
              44},{84,43},{85.875,43}}, color={0,0,127}));
      connect(Alb, combinedModel.Alb) annotation (Line(points={{-1,-67},{18.5,-67},
              {18.5,-66},{40,-66}},color={0,0,127}));
      connect(Pi, combinedModel.Pi) annotation (Line(points={{-1,-83},{18.5,-83},
              {18.5,-74},{40,-74}}, color={0,0,127}));
      connect(ctO2content.pO2, inverseBlockConstraints.y1) annotation (Line(
            points={{38,76},{4,76},{4,62},{-1.4,62}}, color={0,0,127}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={238,46,47},
              pattern=LinePattern.None,
              fillColor={255,255,170},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),
            Rectangle(
              extent={{0,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{0,70},{60,90}},
              lineColor={28,108,200},
              textString="Volume"),
            Text(
              extent={{0,30},{60,50}},
              lineColor={28,108,200},
              textString="T"),
            Text(
              extent={{0,-10},{60,10}},
              lineColor={28,108,200},
              textString="ctHb"),
            Text(
              extent={{0,-50},{60,-30}},
              lineColor={28,108,200},
              textString="Alb"),
            Text(
              extent={{0,-90},{60,-70}},
              lineColor={28,108,200},
              textString="Pi"),
            Polygon(
              points={{-34,40},{-54,-4},{-32,-22},{-12,-6},{-32,40},{-34,40}},
              lineColor={0,0,0},
              fillColor={238,46,47},
              fillPattern=FillPattern.Solid,
              smooth=Smooth.Bezier),
            Ellipse(
              extent={{-44,10},{-20,-14}},
              fillColor={255,255,170},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Ellipse(
              extent={{-48,16},{-20,-12}},
              fillColor={238,46,47},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Text(
              extent={{64,88},{96,70}},
              lineColor={28,108,200},
              fillColor={238,46,47},
              fillPattern=FillPattern.Solid,
              textString="=X0"),
            Text(
              extent={{64,50},{96,32}},
              lineColor={28,108,200},
              fillColor={238,46,47},
              fillPattern=FillPattern.Solid,
              textString="=X0"),
            Text(
              extent={{62,10},{94,-8}},
              lineColor={28,108,200},
              fillColor={238,46,47},
              fillPattern=FillPattern.Solid,
              textString="=X0"),
            Text(
              extent={{62,-30},{94,-48}},
              lineColor={28,108,200},
              fillColor={238,46,47},
              fillPattern=FillPattern.Solid,
              textString="=X0"),
            Text(
              extent={{62,-70},{94,-88}},
              lineColor={28,108,200},
              fillColor={238,46,47},
              fillPattern=FillPattern.Solid,
              textString="=X0"),
            Text(
              extent={{-100,-100},{0,-22}},
              lineColor={238,46,47},
              fillColor={238,46,47},
              fillPattern=FillPattern.Solid,
              textString="FULL
BLOOD"),
            Text(
              extent={{-78,60},{-50,80}},
              lineColor={102,44,145},
              textString="O2"),
            Text(
              extent={{-78,20},{-50,40}},
              lineColor={102,44,145},
              textString="CO2"),
            Text(
              extent={{-78,-20},{-50,0}},
              lineColor={102,44,145},
              textString="BE")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p><span style=\"font-family: Arial,sans-serif; color: #222222; background-color: #ffffff;\">Full-blood model subcomponent with interfaces for complex models. Supplement to the article.</span></p>
</html>"));
    end FullBlood;

    model testFullBloodSubModel
      extends Modelica.Icons.Example;
      FullBlood fullBlood(
        UseVolInput=false,
        UseTInput=false,
        UsectHctInput=false,
        UseAlbInput=false,
        UsePiInput=false,
        Vol0(displayUnit="l"),
        ctHb0=8.4)
        annotation (Placement(transformation(extent={{-20,-40},{60,40}})));
      Physiolibrary.Chemical.Sources.UnlimitedSolutePump unlimitedSolutePump(
          useSoluteFlowInput=false)
        annotation (Placement(transformation(extent={{-60,18},{-40,38}})));
      Physiolibrary.Chemical.Sources.UnlimitedSolutePump unlimitedSolutePump1(
          useSoluteFlowInput=false)
        annotation (Placement(transformation(extent={{-60,2},{-40,22}})));
      Physiolibrary.Chemical.Sources.UnlimitedSolutePump unlimitedSolutePump2(
          useSoluteFlowInput=false, SoluteFlow(displayUnit="mol/s") = 0)
        annotation (Placement(transformation(extent={{-60,-14},{-40,6}})));
    equation
      connect(unlimitedSolutePump.q_out, fullBlood.port_O2) annotation (Line(
          points={{-40,28},{-16,28}},
          color={107,45,134},
          thickness=1));
      connect(fullBlood.port_CO2, unlimitedSolutePump1.q_out) annotation (Line(
          points={{-16,12},{-40,12}},
          color={107,45,134},
          thickness=1));
      connect(fullBlood.port_BE, unlimitedSolutePump2.q_out) annotation (Line(
          points={{-16,-4},{-40,-4}},
          color={107,45,134},
          thickness=1));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p><span style=\"font-family: Arial,sans-serif; color: #222222; background-color: #ffffff;\">Run this model and observe e.g. the pH (fullBlood.combinedModel.combinedModel.pH)</span></p>
<p><br><span style=\"font-family: Arial,sans-serif; color: #222222; background-color: #ffffff;\">This parametrization fits the normal arterial blood. You can run static characteritics by setting the flow into the pumps - either O2, CO2, or HCO3 (in the form of BE).</span></p>
</html>"));
    end testFullBloodSubModel;
  end FullBloodSubmodelComponent;

  package Wolf
    package OriginalValues
      package Auxiliary
        record Erythrocyte
          import Modelica.SIunits.*;

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
          Real Cl_mass0(unit="mol") = Cl0/wf0*water_volume0;
          Real Cl_mass(unit="mol") = Cl/wf0*water_volume;

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
          Concentration water_c[cont]=volume_c ./ wf0 .* water_volume0 ./
              water_volume
            "Actual concentration recalculated to water fraction (initially 0.73) per one liter of volume";
          //   Real masses0[cont](each unit="mol") = volume_c ./ wf0 .* water_volume0
          //     "total mass of contents";

          Real ZHb=15.6 - 23*(10^(pH - 6.69)/(1 + 10^(pH - 6.69))) - 4*(10^(pH
               - 7.89)/(1 + 10^(pH - 7.89)))
            "For spO2=100, otherwise Zhb += 1.5*((1 - 0.75)/0.75)";
          Real ZDPG=(-3) - 1*(10^(pH - 7.56)/(1 + 10^(pH - 7.56))) - 1*(10^(pH
               - 7.32)/(1 + 10^(pH - 7.32)));
          Real ZATP=(-3) - 1*(10^(pH - 6.8)/(1 + 10^(pH - 6.8)));
          Real ZGSH=(-1) - 1*(10^(pH - 8.54)/(1 + 10^(pH - 8.54))) - 1*(10^(pH
               - 9.42)/(1 + 10^(pH - 9.42)));
          Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87))
            "this was not in the original implementatio by  machek!";
          //   Real Zim;//=-9.2 "Charge of ALL impermeable solutes";
          parameter Real Zim=-4.85 "Adjusted Charge of ALL impermeable solutes";

          //osmolality
          Real fiHb=1 + 0.115*water_c[cont.Hb] + 0.0256*water_c[cont.Hb]^2
            "Osmotic coefficient of Hemoglobin (eq from Raftos et al)";
          Real fi[cont]={0.93,0.93,0.93,fiHb,1,1,1,1,0.93}
            "Osmotic coefficients";
          Real Osm=sum(water_c .* fi) + 0.93*(HCO3 + CO3)
            "osmotic coefficient from Wolf model V3.49";
          //   //electric charge
          Real Z[cont]={1,1,-1,ZHb,ZDPG,ZATP,ZGSH,0,ZPi}
            "equivalent charges per mole";
          Real charge=(sum(water_c .* Z) - HCO3 - 2*CO3 + Zim)*water_volume;

          // parameter Fraction Hct0 = 0.44;
          // parameter Volume totalBlodVolume;
          parameter Volume water_volume0=0.44*0.73*1;
          //Hct * wf0 * Vblood;
          constant VolumeFraction wf0=0.73
            "Fraction of water in the erythrocyte";

          Concentration HCO3=0.026*pCO2mmHg*10^(pH - 6.11)
            "HCO3 should be 13.6 acc to article";
          Concentration CO3=HCO3*10^(pH - 10.2);
          Real pH(start=7.22) = -log10(H);

          Volume water_volume(start=water_volume0);
          Real pCO2mmHg(unit="1");
          Real H(start=10^(-7.2));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end Erythrocyte;

        record Plasma
          import Modelica.SIunits.*;

          type cont = enumeration(
              Na,
              K,
              Ca,
              Mg,
              Cl,
              Pi,
              Alb,
              im) "Contents of erythrocyte";

          parameter Concentration Na=140;
          parameter Concentration K=4.1;
          parameter Concentration Ca=2.3;
          parameter Concentration Mg=0.8;
          parameter Concentration Cl0=105;
          Concentration Cl(start=Cl0);
          Real Cl_mass0(unit="mol") = Cl0/wf0*water_volume0;
          Real Cl_mass(unit="mol") = Cl/wf0*water_volume;
          parameter Concentration Pi=1.2;
          parameter Concentration Alb=0.65;
          // parameter Concentration im = 0; // original
          parameter Concentration im=11.87;
          // adjusted

          Concentration volume_c[cont]={Na,K,Ca,Mg,Cl,Pi,Alb,im}
            "concentration in one liter";
          Concentration water_c[cont]=volume_c/wf0*water_volume0/water_volume
            "Actual concentration recalculated to water fraction in one liter (initially 0.73)";
          //   Real masses0[cont](each unit="mol") = volume_c ./ wf0 .* water_volume0
          //     "total mass of contents";

          //charge on inpermeable solutes
          Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87));
          Real ZAlb=(-10.7) - 16*(10^(pH - 7.42)/(1 + 10^(pH - 7.42)));
          //   parameter Real Zim ;//= -5.3 "Charge of ALL impermeable solutes";//test
          // Real Zim = -5.3 "Charge of ALL impermeable solutes";//test
          parameter Real Zim=-10.99
            "ADJUSTED Charge of ALL impermeable solutes";
          //test

          Real fi[cont]={0.93,0.93,1,1,0.93,0.93,1,1};
          // Alb has osmotic coefficient of 1???
          Real Osm=sum(water_c .* fi) + 0.93*(HCO3 + CO3);

          Real Z[cont]={1,1,2,2,-1,ZPi,ZAlb,0};
          Real Chrgs[cont]=water_c .* Z;
          Real charge=(sum(water_c .* Z) + Zim - HCO3 - 2*CO3);
          //*water_volume;

          parameter Volume water_volume0=(1 - 0.44)*0.94*1;

          constant VolumeFraction wf0=0.94 "Fraction of water in the plasma";
          Concentration HCO3=0.0306*pCO2mmHg*10^(pH - 6.11);
          Concentration CO3=HCO3*10^(pH - 10.2);
          Real pH(start=7.37) = -log10(H);

          Volume water_volume(start=water_volume0);
          Real pCO2mmHg(unit="1");
          Real H(start=10^(-7.37));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end Plasma;
      end Auxiliary;

      model PE "Plasma - erythrocyte model"

        FullBloodAcidBase.Wolf.OriginalValues.Auxiliary.Erythrocyte ery;
        FullBloodAcidBase.Wolf.OriginalValues.Auxiliary.Plasma pla;

        // total mass of Cl mobile ion
        Real totalCl=ery.Cl_mass0 + pla.Cl_mass0;
        Real totalWater=ery.water_volume0 + pla.water_volume0;
        Real pCO2mmHg=40;
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
        ery.water_c[ery.cont.Cl]/pla.water_c[pla.cont.Cl] = pla.H/ery.H;

        // pco2
        ery.pCO2mmHg = pCO2mmHg;
        pla.pCO2mmHg = pCO2mmHg;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end PE;

      package Tests
        model P
          FullBloodAcidBase.Wolf.OriginalValues.Auxiliary.Plasma ery;

          // total mass of Cl mobile ion
          Real pCO2mmHg=40;
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

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end P;

        model E
          FullBloodAcidBase.Wolf.OriginalValues.Auxiliary.Erythrocyte ery;

          // total mass of Cl mobile ion
          Real pCO2mmHg=40;
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

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end E;
      end Tests;
    end OriginalValues;

    package v349
      package Auxiliary
        record Erythrocyte
          import Modelica.SIunits.*;

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
          //  Real Cl_mass = Cl/wf0*Vew;

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
          Real HbGperLiterBlood=Hb*Ve0ByVeFraction*fH*64.5;
          Concentration volume_c[cont]={Na,K,Cl,Hb,DPG,ATP,GSH,im,Pi}
            "concentration in one liter";
          Concentration water_c[cont]=volume_c ./ wf0 .* Vew0 ./ Vew
            "Actual concentration recalculated to water fraction (initially 0.73) per one liter of volume";
          //   Real masses0[cont](each unit="mol") = volume_c ./ wf0 .* water_volume0
          //     "total mass of contents";

          Real ZHb=15.6 - 23*(10^(pH - 6.69)/(1 + 10^(pH - 6.69))) - 4*(10^(pH
               - 7.89)/(1 + 10^(pH - 7.89))) + (1 - O2Sat)*HbGperLiterBlood*
              0.012
            "in article is Zhb += 1.5*((1 - 0.75)/0.75) in venous blood";
          Real ZDPG=(-3) - 1*(10^(pH - 7.56)/(1 + 10^(pH - 7.56))) - 1*(10^(pH
               - 7.32)/(1 + 10^(pH - 7.32)));
          Real ZATP=(-3) - 1*(10^(pH - 6.8)/(1 + 10^(pH - 6.8)));
          Real ZGSH=(-1) - 1*(10^(pH - 8.54)/(1 + 10^(pH - 8.54))) - 1*(10^(pH
               - 9.42)/(1 + 10^(pH - 9.42)));
          Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87))
            "this was not in the original implementatio by  machek!";
          //   Real Zim;//=-9.2 "Charge of ALL impermeable solutes";
          //   parameter Real Zim=-4.85 "Adjusted Charge of ALL impermeable solutes";
          parameter Real Zim=-8.0452
            "Charge of ALL impermeable solutes acc to MODEL";

          //osmolality
          Real fiHb=1 + 0.115*water_c[cont.Hb] + 0.0256*water_c[cont.Hb]^2
            "Osmotic coefficient of Hemoglobin (eq from Raftos et al)";
          Real fi[cont]={0.93,0.93,0.93,fiHb,1,1,1,1,1} "Osmotic coefficients";
          Real Osm=sum(water_c .* fi) + 0.93*(HCO3 + CO3) + Lactate +
              permeableParticles
            " = 284.084 osmotic coefficient from Wolf model V3.49";
          Real OsmTest[:]=water_c .* fi;
          //   //electric charge
          Real Z[cont]={1,1,-1,ZHb,ZDPG,ZATP,ZGSH,0,ZPi}
            "equivalent charges per mole";
          Real chargeTest[:]=water_c .* Z;
          Real charge=(sum(water_c .* Z) - HCO3 - 2*CO3 - Lactate + Zim*Vew0/
              Vew)*Vew;

          // parameter Fraction Hct0 = 0.44;
          // parameter Volume totalBlodVolume;
          Real Vew0;
          constant Real wf0(unit="1") = 0.73
            "Fraction of water in the erythrocyte";

          Concentration HCO3=0.026*pCO2mmHg*10^(pH - 6.103)/few
            "HCO3 should be 16,8 acc to the model";
          Concentration CO3=HCO3*10^(pH - 10.2) "CO3 is 0,016";
          Real pH(start=7.22) = -log10(H);

          Real Vew(
            unit="l",
            min=0,
            max=10,
            start=1.43);
          Real pCO2mmHg(unit="1");
          Real H(
            start=10^(-7.2),
            min=0,
            max=1);
          Concentration Lactate "LACew =LACpw/rCl";
          parameter Concentration permeableParticles=10.64
            "glucose and urea concentration in PLasma water";

          Concentration HCO3ClDiff=HCO3/water_c[cont.Cl]*few;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end Erythrocyte;

        record Plasma
          import Modelica.SIunits.*;

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

          Concentration Ca=2.4;
          parameter Concentration Mg=0.8;
          //parameter Concentration Cl0=105;
          //  Real Cl_mass0(unit="mol") = Cl0/wf0*water_volume0;
          //  Real Cl_mass(unit="mol") = Cl/wf0*water_volume;
          parameter Concentration Pi=1.2;
          Concentration Alb=AlbPwGPerL/66.5 "mmol/Lpw";
          Real AlbPwGPerL=43*Vp0ByVp "g/lpw";
          // parameter Concentration im = 0; // original
          parameter Concentration im=11.87;
          // adjusted
          parameter Concentration Lac=1.5;
          Concentration SO4pw=0.33/fpw "mmol/Lpw";
          Concentration volume_c[cont]={Na,K,Ca,Mg,Cl,Pi,Alb/Vp0ByVp,im,Lac}
            "concentration in one liter";
          Concentration water_c[cont]=volume_c/fpw*Vp0ByVp
            "Actual concentration recalculated to water fraction in one liter (initially 0.96)";
          //   Real masses0[cont](each unit="mol") = volume_c ./ wf0 .* water_volume0
          //     "total mass of contents";
          Real CaPPBound=(3.98)/(3.98 + Hh)*water_c[cont.Ca] "1.257";
          Real MgPPBound=CaPPBound/2;
          Real ZCaBindPerAlb=CaPPBound/water_c[cont.Alb];
          Real ZClBindPerAlb=6.46/(6.46 + Hh)*6 + 4 "Anstey";
          Real CaIon=water_c[cont.Ca]*2 - CaPPBound "mEq/Lpw";
          Real MgIon=water_c[cont.Mg]*2 - MgPPBound "mEq/Lpw";
          //charge on inpermeable solutes
          Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87));
          Real ZFigge=(-10.65) - 16*(10^(pH - 7.418)/(1 + 10^(pH - 7.418)));
          Real ZAlbBnd=-ZClBindPerAlb + ZCaBindPerAlb + ZCaBindPerAlb/2;
          Real ZAlb=ZFigge + ZAlbBnd;
          Real Hh=H/fpw*1e8;

          //   parameter Real Zim ;//= -5.3 "Charge of ALL impermeable solutes";//test
          // Real Zim = -5.3 "Charge of ALL impermeable solutes";//test
          parameter Real Zim=-0.6349 "Charge of ALL impermeable solutes";

          Real fi[cont]={0.93,0.93,1,1,0.93,0.93,1,1,-9999};
          // Alb has osmotic coefficient of 1???
          Real CaOsm=water_c[cont.Ca] - CaPPBound*0.5;
          Real MgOsm=water_c[cont.Mg] - MgPPBound*0.5;
          Real OsmPart=SO4pw + CaOsm + MgOsm + sum(water_c[{cont.Na,cont.K,cont.Cl,
              cont.Lac}]) + water_c[cont.Pi] + HCO3 + CO3;
          Real Osm=OsmPart*0.93 + permeableParticles " = 284.084";

          Real Z[cont]={1,1,2,2,-1,ZPi,ZAlb,0,1};
          Real Chrgs[cont]=water_c .* Z;

          Real SID=sum(water_c[{cont.Na,cont.K}]) - water_c[cont.Cl] + 2*sum(
              water_c[{cont.Ca,cont.Mg}]) - CaPPBound - MgPPBound - water_c[
              cont.Lac];
          Real charge=Zim + SID + water_c[cont.Alb]*ZAlb + water_c[cont.Pi]*ZPi
               - HCO3 - 2*CO3 - SO4pw*2;
          //*water_volume;

          parameter Volume water_volume0=(1 - 0.44)*0.94*1;

          constant VolumeFraction wf0=0.94 "Fraction of water in the plasma";
          Concentration HCO3=0.0306*pCO2mmHg*10^(pH - 6.103)/fpw;
          Concentration CO3=HCO3*10^(pH - 10.2);
          Real pH(start=7.37) = -log10(H);

          Real Vp0ByVp;
          Real fpw "fraction plasma water - Vpw/Vp";
          Real pCO2mmHg(unit="1");
          Real H(
            start=10^(-7.37),
            min=0,
            max=1);
          Real Hw=H/fpw;
          Concentration permeableParticles=(5 + 5)/fpw
            "glucose and urea concentration in PLasma water";
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end Plasma;

        record Isf
          "ISF is considered as water compartment, that is each substance concentration is meant as in water"
          import Modelica.SIunits.*;

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
          Concentration MgIon "mEq/l";
          //Modelica.SIunits.Concentration CaPPBound "mEq/l";
          Concentration Cl;
          Real Vis;
          Real Vis0;
          parameter Real Alb(unit="g/l") = 43 "plasma alb conc g/l";
          Concentration SO4pw "mmol/Pw";

          Concentration albis=(Alb*0.33*Vis0)/(Vis - Vis0*0.25) "g/iw";
          Concentration albiw=albis/66.5
            "mmol/liw Alb conc is fixed as 1/3 of 43, no matter the conc in plasma";

          Real CaPPBound=(3.98)/(3.98 + Hh)*plasma_water_c[cont.Ca]
            "= 1.290 Meq/Lis";
          Real CaPPBoundErr=0
            "Error in VisSim diagram, should be CaPPBound instead";
          Real ZCaBindPerAlb=CaPPBound/albiw;
          Real ZClBindPerAlb=6.46/(6.46 + Hh)*6 + 4 "Anstey";

          Real transf[cont]={rClpwis,rClpwis,0,0,0,(1/rClpwis)^(-ZPi),0,0,1/
              rClpwis} "Only Na, K, Pi and Lac";

          Real SO4=(SO4pw*2)/(rClpwis^2) "concentration in one liter mEq/lw";
          // Concentration water_c[cont]=transf .* plasma_water_c;
          Concentration water_c[cont]=transf .* plasma_water_c
            "Actual concentration recalculated from plasma";

          Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87));
          Real ZFigge=(-10.65) - 16*(10^(pH - 7.418)/(1 + 10^(pH - 7.418)));
          Real ZAlbBnd=-ZClBindPerAlb + ZCaBindPerAlb + ZCaBindPerAlb/2;
          Real ZAlb=ZFigge + ZAlbBnd;

          parameter Real Zim=3.057 "Charge of ALL impermeable solutes";

          Real OsmStep1=sum(water_c[{cont.Na,cont.K}]) + Cl;
          Real CaMgIonOsm=(CaIon - CaPPBoundErr + MgIon - CaPPBoundErr/2)/2;
          Real OsmStep2=OsmStep1 + CaMgIonOsm + sum(water_c[{cont.Pi,cont.Lac}])
               + HCO3 + CO3 + SO4/2;
          Real Osm=OsmStep2*0.93 + permeableParticles;

          Real SID=sum(water_c[{cont.Na,cont.K}]) - Cl + CaIon - CaPPBoundErr
               + (MgIon - 0.5*CaPPBoundErr) - water_c[cont.Lac];
          Real charge=Zim + SID + albiw*ZAlb + water_c[cont.Pi]*ZPi - HCO3 - 2*
              CO3 - SO4;
          //*water_volume;

          Concentration HCO3=0.0326*pCO2mmHg*10^(pH - 6.103);
          Concentration CO3=HCO3*10^(pH - 10.2);
          Real pH(start=7.408) = -log10(H);

          Real pCO2mmHg(unit="1");
          Real H(
            start=10^(-7.408),
            min=0,
            max=1);
          Real Hh=H*1e8;

          parameter Concentration permeableParticles=10.64
            "glucose and urea concentration in PLasma water";

          // Transcapillary pressures
          // wattenpaugh J trauma Inuj 1998
          Real AlbPwGPerL;
          Real ca=AlbPwGPerL*0.1 "Albumin concentration";
          Real cp=AlbPwGPerL*0.1*1.637 "protein concentration";
          Real COPPpl=PIEalb + PIEglob;
          Real PIEalb=ca/cp*(cp*2.8 + 0.18*cp^2 + 0.012*cp^3);
          Real PIEglob=(1 - ca/cp)*(0.9*cp + 0.12*cp^2 + 0.004*cp^3);
          Real cai=albis*0.1;
          Real COPalbis=cai*2.8 + 0.18*cai^2 + 0.012*cai^3;
          Real BloodyHell=(Vis - Vis0)/Vis0;
          Real Bh1=if BloodyHell < 0 then BloodyHell*125 else BloodyHell*243
            "Simplified bloody hell";
          //   Real pd = (28.48 - 0.99*(COPPpl - COPalbis))/19.3;
          Real pd=(38.36 - 0.99*(COPPpl - COPalbis) + PVb - Bh1)/19.3;
          Real OsmDiff=pd*2;
          Real Vbr=(Vb - Vb0)/Vb0 "VbRatio";
          Real PVb=if Vbr < -5 then -5*50 elseif Vbr > -5 and Vbr < 0 then Vbr*
              50 elseif Vbr > 0 and Vbr < 100 then Vbr*46.5 else 100*46.5
            "In original -5..0 and 0..100";
          Real Vb;
          Real Vb0;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end Isf;

        record Volumes
          type PatientCat = enumeration(
              adultM,
              adultF,
              childM,
              childF);
          parameter PatientCat pc=PatientCat.adultM;
          parameter Real Fcell=if pc == PatientCat.adultM or pc == PatientCat.adultF
               then 0.88 elseif pc == PatientCat.childM then 0.83 else 0.85;

          // the rest is only adult male
          parameter Real TBWm=2.447 - 0.09156*age + 0.3362*weight + 0.1074*
              height;
          //Real TBWf
          parameter Real LBM=0.407*weight + 0.267*height - 19.2;
          // Real LBMm, LBMf;
          parameter Real age=21;
          parameter Real weight=70;
          parameter Real height=175;

          parameter Real Ve0=0.048*LBM - 0.72;
          parameter Real Vp0=0.051*LBM + 0.27;
          parameter Real Vb0=(Ve0 + Vp0);
          parameter Real fH0=Ve0/(Ve0 + Vp0);
          parameter Real fHa0=fH0/Fcell;
          parameter Real Vis0=0.331*LBM - Vpw0;
          parameter Real TBW=2.447 - 0.09156*age + 0.3362*weight + 0.1074*
              height "Woodrow 2003 for male, there also exists female variant";
          parameter Real Vc0=TBW - Vis0 - Vpw0 - Vew0;

          parameter Real Vpw0=Vp0 - Vps;
          parameter Real Vps=0.0603*Vp0;
          parameter Real Vew0=Ve0*0.7317;
          parameter Real Ves=Ve0 - Vew0;
          Real Ve=Ves + Vew;
          Real Vb=Vpl + Ve;
          Real fH=Ve/Vb;

          parameter Real Vw0=Ve0 + Vp0 + Vis0 + Vc0;
          parameter Real V=Vw0 + Vadd;
          parameter Real Vadd=0;
          // Volumes
          Real Vpl=Vpw + Vps;
          Real fpw=Vpw/Vpl;
          Real few=Vew/Ve;
          Real Vpw=Vadd + Vpw0 + Vis0 + Vew0 + Vc0 - Vis - Vew - Vc " = 2.93";
          Real Vp0ByVp=Vp0/Vpl;
          Real Vew(start=Vew0,min=1e-3) "= 1.43";
          Real Vis(start=Vis0,min=1e-3) "= 15.59";
          Real Vc(start=Vc0,min=1e-3) "= 22.9";
          // Real Vew = 1.43;
          // Real Vis = 15.59;
          // Real Vc = 22.9;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end Volumes;

        record StrongIonMasses

          Real Ve0;
          Real Vp0;
          Real Vis0;
          Real Vc0;
          // any ion addition MUST be added here as well;
          constant Real Napl0=140;
          constant Real Nais0=141.5;
          constant Real Nac0=12;
          Real MNa=Napl0*Vp0 + Nais0*Vis0 + Nac0*Vc0;
          Real MNac=Nac0*Vc0;
          Real MNa_IP;

          constant Real Kpl0=4.7;
          constant Real Kis0=4.75;
          constant Real Kc0=139;
          Real MK=Kpl0*Vp0 + Kis0*Vis0 + Kc0*Vc0;
          Real MKc=Kc0*Vc0;
          Real MK_IP;

          constant Real Cle0=53.5;
          constant Real Clpl0=104;
          constant Real Clis0=116.5;
          constant Real Clc0=4;
          Real MCl=Cle0*Ve0 + Clpl0*Vp0 + Clis0*Vis0 + Clc0*Vc0;
          Real MClc=Clc0*Vc0;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end StrongIonMasses;

        record Cell
          import Modelica.SIunits.*;

          Real Kis;
          Real Nais;
          Real Clis;
          Real Vc0;
          Real Vc;
          Real His;
          Real pCO2mmHg(unit="1");

          Real Clc(
            start=3.9,
            min=1,
            max=100);

          Real rClisce(start=29.6) = Clis/Clc;
          Real Kc=Kis*rClisce;
          Real Nac=Nais*rClisce*0.0029;

          Concentration HCO3=0.029*pCO2mmHg*10^(pH - 6.103);
          Real pH(start=7.408) = -log10(H);

          Real H(
            start=10^(-6.935),
            min=1e-9,
            max=1) = His*0.1*rClisce "Why His*0.1? Is it the water content?";

          Real charge=-HCO3 + Kc - Clc + Nac + Zim*Im;
          Real Zim=-1.1728*10^(pH - 5.5)/(1 + 10^(pH - 5.5));
          Real Im=Vc0/Vc*123.045;

          Real OsmPart=HCO3 + Kc + Nac + Clc;
          Real Osm=Im + permeableParticles + 0.93*OsmPart;

          Real MNac=Nac*Vc;
          Real MClc=Clc*Vc;
          Real MKc=Kc*Vc;

          parameter Concentration permeableParticles=5.313 + UreaError
            "glucose and urea concentration in PLasma water";
          constant Real UreaError=0
            "Concentration of Urea in plasma water. Should be 5.31, is zero";

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end Cell;
      end Auxiliary;

      package Tests
        model E
          FullBloodAcidBase.Wolf.v349.Auxiliary.Erythrocyte ery;
          FullBloodAcidBase.Wolf.v349.Auxiliary.Volumes vols;
          // total mass of Cl mobile ion
          Real pCO2mmHg=40;
          // output: osmolarity and Cl
          parameter Real Hct=0.44;
          parameter Real O2s=0.75;

          parameter Real Hpw=4.128159e-8;
          parameter Real Clpw=109.694599;
          parameter Real EryOsm=284.084;
          Real Cle=ery.water_c[ery.cont.Cl] "= 71.031139 @v3.49";
          Real Vew=vols.Vew "= 1.44400";
        equation
          //  vols.Vew = 1.44400;
          vols.Vis = 15.75516;
          vols.Vc = 22.64963;

          vols.Ve0/vols.Ve = ery.Ve0ByVeFraction;
          ery.fH = vols.fH;
          ery.few = vols.few;
          ery.Lactate = 1.0185;
          //LACpw/rCl;
          ery.Vew0 = vols.Vew0;
          ery.Vew = vols.Vew;
          Cle/Clpw = Hpw/ery.H "Clpla , Hpl at standard V3.49";

          // ery.water_c[ery.cont.Cl] = 71.0311;
          // ery.H = 6.37518e-8;
          // ery.
          //   // same osmolarities
          ery.Osm = EryOsm;
          ery.charge = 0;
          //   ery.Cl = ery.Cl0;
          //ery.Cl = 71.77 * ery.wf0;

          // Donnan equilibrium
          // pco2
          ery.pCO2mmHg = pCO2mmHg;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end E;

        model P
          Real Clpw(
            start=112,
            min=10,
            max=200) = pla.water_c[pla.cont.Cl] " = 109.694599";
          Real Vis=vols.Vis "=15.75516";

          FullBloodAcidBase.Wolf.v349.Auxiliary.Plasma pla;
          FullBloodAcidBase.Wolf.v349.Auxiliary.Volumes vols;
          FullBloodAcidBase.Wolf.v349.Auxiliary.StrongIonMasses sim;
          // total mass of Cl mobile ion
          Real pCO2mmHg=40;
          // Real rClpw_is = pla.volume_c[pla.cont.Cl]/Clis;
          Real rClpwis=Clpw/Clis;

          constant Real Clis=115.461419;
          constant Real MNac=272.17930;
          constant Real MKc=3178.3177;
          constant Real MClc=88.33349;
          constant Real MCle=102.5695;

          parameter Real EryOsm=284.084;
          Real MCl_IP=sim.MCl - MClc - MCle;
          //Real test = sim.MNa_IP /(vols.Vis * rClpw_is + vols.Vpw);
        equation
          vols.Vew = 1.44400;
          //  vols.Vis = 15.75516;
          vols.Vc = 22.64963;

          sim.Ve0 = vols.Ve0;
          sim.Vp0 = vols.Vp0;
          sim.Vis0 = vols.Vis0;
          sim.Vc0 = vols.Vc0;
          sim.MNa_IP = sim.MNa - MNac;
          sim.MK_IP = sim.MK - MKc;

          pla.water_c[pla.cont.Na] = sim.MNa_IP/(vols.Vis*rClpwis + vols.Vpw);
          pla.water_c[pla.cont.K] = sim.MK_IP/(vols.Vis*rClpwis + vols.Vpw);
          Clis = (MCl_IP - vols.Vpw*pla.water_c[pla.cont.Cl])/vols.Vis;
          pla.water_c[pla.cont.Cl] = 109.694599;

          // same osmolarities
          //  pla.Osm = EryOsm;
          pla.charge = 0;

          pla.Vp0ByVp = vols.Vp0ByVp;

          // Donnan equilibrium
          //  pla.pH = 7.410676;
          pla.fpw = vols.fpw;
          // pco2
          pla.pCO2mmHg = pCO2mmHg;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end P;

        model I
          FullBloodAcidBase.Wolf.v349.Auxiliary.Isf isf(permeableParticles=
                10.62755);
          FullBloodAcidBase.Wolf.v349.Auxiliary.Volumes vols;
          Real rClpwis=Clpw/isf.Cl " = 0.9500";
          parameter Real Hpw=4.128159e-8;
          parameter Real Clpw=109.694599;
          parameter Real Alb=43;
          constant Real EryOsm=284.084;
          constant Real MCl_IP=2148.68379;

          // Real Clis( start = 115) = MCl_IP - vols.Vpw*Clpw/vols.Vis " = 115.4614";
          Real Clis(start=115) " = 115.4614";
          Real Vis=vols.Vis "= 15.75516";
          Real oo=isf.Osm + isf.pd*2;
        equation
          vols.Vew = 1.44400;
          //  vols.Vis = 15.75516;
          vols.Vc = 22.64963;

          isf.plasma_water_c = {147.3244,4.9890,2.4977,0,0,1.2488,0,0,1.56017};

          isf.rClpwis = rClpwis;

          isf.CaIon = 3.76940*rClpwis^2;
          isf.MgIon = 1.05212*rClpwis^2;
          isf.Vis = vols.Vis;
          isf.Vis0 = vols.Vis0;
          isf.SO4pw = 0.702420/2;
          // The input here is in mmol/lpw to match plaswma variable, not in mEq/lpw as in the vissim model
          isf.pCO2mmHg = 40;
          isf.AlbPwGPerL = Alb*vols.Vp0ByVp;
          isf.Vb = vols.Vb;
          isf.Vb0 = vols.Vb0;
          isf.Cl = Clis;

          //  isf.permeableParticles = 10.62755 "5+5 / fpw";

          //  isf.Cl = 115.4614;
          //  isf.pH = 7.40649;

          isf.H = Hpw*rClpwis;
          isf.charge = 0;
          isf.Osm + isf.OsmDiff = EryOsm;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end I;

        model C
          import Modelica.SIunits.*;
          FullBloodAcidBase.Wolf.v349.Auxiliary.Cell c;
          Real Clc(
            start=3.9,
            min=1,
            max=100) = c.Clc;
        equation
          c.Kis = 4.73983;
          //4.73106;
          c.Nais = 139.9662;
          //141.474;
          c.Clis = 115.461419;
          //116.7911;
          c.Vc0 = 22.87176;
          //  c.Vc = 22.64963;//22.897789;
          c.His = 3.92197e-8;
          //3.91020e-8;
          c.pCO2mmHg = 40;

          //  c.rClisce = 29.6055;
          //  c.Clc = 3.961488;
          c.charge = 0;
          c.Osm = 282.20;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end C;
      end Tests;

      model EP

        Real Clpw(
          start=112,
          min=10,
          max=200) = pla.water_c[pla.cont.Cl] " = 109.694599";
        Real Cle(
          start=71,
          min=1,
          max=200) = ery.water_c[ery.cont.Cl] "= 71.031139 @v3.49";
        Real Vew=vols.Vew "= 1.44400";
        Real Vis=vols.Vis "=15.75516";

        Real pHe=ery.pH;
        Real pHp=pla.pH;

        FullBloodAcidBase.Wolf.v349.Auxiliary.Plasma pla;
        FullBloodAcidBase.Wolf.v349.Auxiliary.StrongIonMasses sim;
        FullBloodAcidBase.Wolf.v349.Auxiliary.Erythrocyte ery;
        FullBloodAcidBase.Wolf.v349.Auxiliary.Volumes vols;
        // total mass of Cl mobile ion
        Real pCO2mmHg=40;
        // output: osmolarity and Cl
        parameter Real Hct=0.44;
        parameter Real O2s=0.75;

        Real rClpwis=Clpw/Clis;

        constant Real Clis=115.461419;
        constant Real MNac=272.17930;
        constant Real MKc=3178.3177;
        constant Real MClc=88.33349;
        Real MCle=ery.water_c[ery.cont.Cl]*vols.Vew "= 102.5695";
        Real MCl_IP=sim.MCl - MClc - MCle;

      equation
        //  vols.Vew = 1.44400;
        //  vols.Vis = 15.75516;
        vols.Vc = 22.64963;

        vols.Ve0/vols.Ve = ery.Ve0ByVeFraction;
        ery.fH = vols.fH;
        ery.few = vols.few;
        ery.Lactate = pla.water_c[pla.cont.Lac]*Cle/Clpw;
        //LACpw/rCl;
        ery.Vew0 = vols.Vew0;
        ery.Vew = vols.Vew;

        Cle/Clpw = pla.Hw/ery.H "Clpla , Hpl at standard V3.49";
        ery.Osm = pla.Osm;
        ery.charge = 0;
        ery.pCO2mmHg = pCO2mmHg;

        sim.Ve0 = vols.Ve0;
        sim.Vp0 = vols.Vp0;
        sim.Vis0 = vols.Vis0;
        sim.Vc0 = vols.Vc0;
        sim.MNa_IP = sim.MNa - MNac;
        sim.MK_IP = sim.MK - MKc;

        pla.water_c[pla.cont.Na] = sim.MNa_IP/(vols.Vis*rClpwis + vols.Vpw);
        pla.water_c[pla.cont.K] = sim.MK_IP/(vols.Vis*rClpwis + vols.Vpw);
        Clis = (MCl_IP - vols.Vpw*pla.water_c[pla.cont.Cl])/vols.Vis;
        pla.water_c[pla.cont.Cl] = 109.694599;

        pla.charge = 0;
        pla.Vp0ByVp = vols.Vp0ByVp;
        pla.fpw = vols.fpw;
        pla.pCO2mmHg = pCO2mmHg;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end EP;

      model EPI

        Real Clis(
          start=115,
          min=10,
          max=200) = isf.Cl " = 115.4614";
        Real Clpw(
          start=112,
          min=10,
          max=200) = pla.water_c[pla.cont.Cl] " = 109.694599";
        Real Cle(
          start=71,
          min=1,
          max=200) = ery.water_c[ery.cont.Cl] "= 71.031139 @v3.49";

        Real Vew=vols.Vew "= 1.44400";
        Real Vis=vols.Vis "=15.75516";

        Real pHe=ery.pH;
        Real pHp=pla.pH;
        Real pHi=isf.pH;

        FullBloodAcidBase.Wolf.v349.Auxiliary.Plasma pla;
        FullBloodAcidBase.Wolf.v349.Auxiliary.StrongIonMasses sim;
        FullBloodAcidBase.Wolf.v349.Auxiliary.Erythrocyte ery;
        FullBloodAcidBase.Wolf.v349.Auxiliary.Isf isf(permeableParticles=
              10.62755);
        FullBloodAcidBase.Wolf.v349.Auxiliary.Volumes vols;
        // total mass of Cl mobile ion
        Real pCO2mmHg=40;
        // output: osmolarity and Cl
        parameter Real Hct=0.44;
        parameter Real O2s=0.75;

        constant Real MNac=272.17930;
        constant Real MKc=3178.3177;
        constant Real MClc=88.33349;
        Real MCle=ery.water_c[ery.cont.Cl]*vols.Vew "= 102.5695";
        Real MCl_IP=sim.MCl - MClc - MCle;
        Real rH=isf.H/pla.Hw;

        Real rClpwis=Clpw/Clis " = 0.9500";
        parameter Real Alb=43;

        // Real Clis( start = 115) = MCl_IP - vols.Vpw*Clpw/vols.Vis " = 115.4614";
      equation
        //  vols.Vew = 1.44400;
        //  vols.Vis = 15.75516;
        vols.Vc = 22.64963;

        vols.Ve0/vols.Ve = ery.Ve0ByVeFraction;
        ery.fH = vols.fH;
        ery.few = vols.few;
        ery.Lactate = pla.water_c[pla.cont.Lac]*Cle/Clpw;
        //LACpw/rCl;
        ery.Vew0 = vols.Vew0;
        ery.Vew = vols.Vew;

        Cle/Clpw = pla.Hw/ery.H "Clpla , Hpl at standard V3.49";
        ery.Osm = pla.Osm;
        ery.charge = 0;
        ery.pCO2mmHg = pCO2mmHg;

        sim.Ve0 = vols.Ve0;
        sim.Vp0 = vols.Vp0;
        sim.Vis0 = vols.Vis0;
        sim.Vc0 = vols.Vc0;
        sim.MNa_IP = sim.MNa - MNac;
        sim.MK_IP = sim.MK - MKc;

        pla.water_c[pla.cont.Na] = sim.MNa_IP/(vols.Vis*rClpwis + vols.Vpw);
        pla.water_c[pla.cont.K] = sim.MK_IP/(vols.Vis*rClpwis + vols.Vpw);
        Clis = (MCl_IP - vols.Vpw*pla.water_c[pla.cont.Cl])/vols.Vis;
        //pla.water_c[pla.cont.Cl] = 109.694599;

        pla.charge = 0;
        pla.Vp0ByVp = vols.Vp0ByVp;
        pla.fpw = vols.fpw;
        pla.pCO2mmHg = pCO2mmHg;

        isf.plasma_water_c = {147.3244,4.9890,2.4977,0,0,1.2488,0,0,1.56017};
        isf.rClpwis = rClpwis;
        isf.CaIon = 3.76940*rClpwis^2;
        isf.MgIon = 1.05212*rClpwis^2;
        isf.Vis = vols.Vis;
        isf.Vis0 = vols.Vis0;
        isf.SO4pw = pla.SO4pw;
        // The input here is in mmol/lpw to match plaswma variable, not in mEq/lpw as in the vissim model
        isf.pCO2mmHg = 40;
        isf.AlbPwGPerL = Alb*vols.Vp0ByVp;
        isf.Vb = vols.Vb;
        isf.Vb0 = vols.Vb0;

        isf.H = pla.Hw*rClpwis;
        isf.charge = 0;
        isf.Osm + isf.OsmDiff = pla.Osm;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end EPI;

      model EPIC
        Real Clc(
          start=3.9,
          min=1,
          max=100) = c.Clc;
        Real Clis(
          start=115,
          min=10,
          max=200) = isf.Cl " = 115.4614";
        Real Clpw(
          start=112,
          min=10,
          max=200) = pla.water_c[pla.cont.Cl] " = 109.694599";
        Real Cle(
          start=71,
          min=1,
          max=200) = ery.water_c[ery.cont.Cl] "= 71.031139 @v3.49";

        Real Vew=vols.Vew "= 1.44400";
        Real Vis=vols.Vis "=15.75516";

        Real pHe=ery.pH;
        Real pHp=pla.pH;
        Real pHi=isf.pH;
        Real pHc=c.pH;

        FullBloodAcidBase.Wolf.v349.Auxiliary.Plasma pla;
        FullBloodAcidBase.Wolf.v349.Auxiliary.StrongIonMasses sim;
        FullBloodAcidBase.Wolf.v349.Auxiliary.Erythrocyte ery;
        FullBloodAcidBase.Wolf.v349.Auxiliary.Isf isf(permeableParticles=
              10.62755);
        FullBloodAcidBase.Wolf.v349.Auxiliary.Volumes vols;
        FullBloodAcidBase.Wolf.v349.Auxiliary.Cell c;
        // total mass of Cl mobile ion
        Real pCO2mmHg=40;
        // output: osmolarity and Cl
        parameter Real Hct=0.44;
        parameter Real O2s=0.75;

        //   constant Real MNac=272.17930;
        //   constant Real MKc=3178.3177;
        //   constant Real MClc=88.33349;
        //   Real MNac=272.17930;
        //   Real MKc=3178.3177;
        //   Real MClc=88.33349;
        Real MNac=c.Nac*vols.Vc;
        Real MKc=c.Kc*vols.Vc;
        Real MClc=c.Clc*vols.Vc;
        Real MCle=ery.water_c[ery.cont.Cl]*vols.Vew "= 102.5695";
        Real MCl_IP=sim.MCl - MClc - MCle;
        Real rH=isf.H/pla.Hw;

        Real rClpwis=Clpw/Clis " = 0.9500";
        parameter Real Alb=43;

        // Real Clis( start = 115) = MCl_IP - vols.Vpw*Clpw/vols.Vis " = 115.4614";
      equation
        //  vols.Vew = 1.44400;
        //  vols.Vis = 15.75516;
        //  vols.Vc = 22.64963;

        vols.Ve0/vols.Ve = ery.Ve0ByVeFraction;
        ery.fH = vols.fH;
        ery.few = vols.few;
        ery.Lactate = pla.water_c[pla.cont.Lac]*Cle/Clpw;
        //LACpw/rCl;
        ery.Vew0 = vols.Vew0;
        ery.Vew = vols.Vew;

        Cle/Clpw = pla.Hw/ery.H "Clpla , Hpl at standard V3.49";
        ery.Osm = pla.Osm;
        ery.charge = 0;
        ery.pCO2mmHg = pCO2mmHg;

        sim.Ve0 = vols.Ve0;
        sim.Vp0 = vols.Vp0;
        sim.Vis0 = vols.Vis0;
        sim.Vc0 = vols.Vc0;
        sim.MNa_IP = sim.MNa - MNac;
        sim.MK_IP = sim.MK - MKc;

        pla.water_c[pla.cont.Na] = sim.MNa_IP/(vols.Vis*rClpwis + vols.Vpw);
        pla.water_c[pla.cont.K] = sim.MK_IP/(vols.Vis*rClpwis + vols.Vpw);
        Clis = (MCl_IP - vols.Vpw*pla.water_c[pla.cont.Cl])/vols.Vis;
        //pla.water_c[pla.cont.Cl] = 109.694599;

        pla.charge = 0;
        pla.Vp0ByVp = vols.Vp0ByVp;
        pla.fpw = vols.fpw;
        pla.pCO2mmHg = pCO2mmHg;

        isf.plasma_water_c = {147.3244,4.9890,2.4977,0,0,1.2488,0,0,1.56017};
        isf.rClpwis = rClpwis;
        isf.CaIon = 3.76940*rClpwis^2;
        isf.MgIon = 1.05212*rClpwis^2;
        isf.Vis = vols.Vis;
        isf.Vis0 = vols.Vis0;
        isf.SO4pw = pla.SO4pw;
        // The input here is in mmol/lpw to match plaswma variable, not in mEq/lpw as in the vissim model
        isf.pCO2mmHg = 40;
        isf.AlbPwGPerL = Alb*vols.Vp0ByVp;
        isf.Vb = vols.Vb;
        isf.Vb0 = vols.Vb0;

        isf.H = pla.Hw*rClpwis;
        isf.charge = 0;
        isf.Osm + isf.OsmDiff = pla.Osm;

        //    c.Kis = 4.73983;//4.73106;
        //    c.Nais = 139.9662;//141.474;
        //    c.Clis = 115.461419;//116.7911;
        //   c.Vc0 = 22.87176;
        //  c.Vc = 22.64963;//22.897789;
        //   c.His = 3.92197e-8;//3.91020e-8;
        //   c.pCO2mmHg = 40;

        //   c.Vc = vols.Vc;
        //  c.rClisce = 29.6055;
        //  c.Clc = 3.961488;
        //   c.charge = 0;
        //   c.Osm = 282.20;
        //

        c.Kis = isf.water_c[isf.cont.K];
        //4.73983;//4.73106;
        c.Nais = isf.water_c[isf.cont.Na];
        //139.9662;//141.474;
        c.Clis = isf.Cl;
        //115.461419;//116.7911;
        c.Vc0 = vols.Vc0;
        c.Vc = vols.Vc;
        c.His = isf.H;
        //3.92197e-8;//3.91020e-8;
        c.pCO2mmHg = 40;

        c.charge = 0;
        c.Osm = isf.Osm;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end EPIC;

      annotation (Documentation(info="<html>
<p>References:</p>
<p><br>Raftos, J. E., Bulliman, B. T., &amp; Kuchel, P. W. (1990). Evaluation of an electrochemical model of erythrocyte pH buffering using 31P nuclear magnetic resonance data. The Journal of General Physiology, 95(6), 1183&ndash;1204.</p>
<p>Wolf, M. B. (2013). Whole body acid-base and fluid-electrolyte balance: a mathematical model. American Journal of Physiology. Renal Physiology, 305(8), F1118&ndash;31.</p>
<p>Wolf, M. B. (2015). Comprehensive diagnosis of whole-body acid-base and fluid-electrolyte disorders using a mathematical model and whole-body base excess. Journal of Clinical Monitoring and Computing, 29(4), 475&ndash;490.</p>
<p>Wolf, M. B., &amp; DeLand, E. C. (2011). A comprehensive, computer-model-based approach for diagnosis and treatment of complex acid&ndash;base disorders in critically-ill patients. Journal of Clinical Monitoring and Computing, 25(6), 353&ndash;364.</p>
<p>Wolf, M. B., &amp; Deland, E. C. (2011). A mathematical model of blood-interstitial acid-base balance: application to dilution acidosis and acid-base status. Journal of Applied Physiology, 110(4), 988&ndash;1002.</p>
</html>"));
    end v349;

    package v351
      package Auxiliary
        record Erythrocyte
          import Modelica.SIunits.*;

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
          parameter Concentration K=130*wf0;
          parameter Concentration Cl0=53.8;
          Concentration Cl(start=Cl0);
          //Real Cl_mass0(unit="mol") = Cl0/wf0*Vew0;
          //  Real Cl_mass = Cl/wf0*Vew;

          Concentration Hb=5.3/(1/wf0*Vew0/Vew)/few
            "concentration of Hb tetramer";
          Concentration DPG=4.3/few*wf0 "the 4.3 value goes directly in";
          Concentration ATP=1.8/few*wf0;

          Concentration GSH=2.23/few*wf0;
          parameter Concentration im=0;
          //16.03*wf0 "original value in model";
          //   parameter Concentration im=21.64 "adjusted value";
          // Concentration im=20.2 "Original value in article [1]";
          //   Concentration im;//=20.2;
          parameter Concentration Pi=0.67 "mono- and di- valent";
          parameter Real O2Sat(unit="1") = 0.75;
          Real Ve0ByVeFraction;
          Real fH;
          Real few;
          Real HbGperLiterBlood=Hb*Ve0ByVeFraction*fH*64.5;
          Concentration volume_c[cont]={Na,K,Cl,Hb,DPG,ATP,GSH,im,Pi}
            "concentration in one liter";
          Concentration water_c[cont]=volume_c ./ wf0 .* Vew0 ./ Vew
            "Actual concentration recalculated to water fraction (initially 0.73) per one liter of volume";
          //   Real masses0[cont](each unit="mol") = volume_c ./ wf0 .* water_volume0
          //     "total mass of contents";

          Real ZHb=15.6 - 23*(10^(pH - 6.69)/(1 + 10^(pH - 6.69))) - 4*(10^(pH
               - 7.89)/(1 + 10^(pH - 7.89))) + (1 - O2Sat)*HbGperLiterBlood*
              0.012
            "in article is Zhb += 1.5*((1 - 0.75)/0.75) in venous blood";
          Real ZDPG=(-3) - 1*(10^(pH - 7.56)/(1 + 10^(pH - 7.56))) - 1*(10^(pH
               - 7.32)/(1 + 10^(pH - 7.32)));
          Real ZATP=(-3) - 1*(10^(pH - 6.8)/(1 + 10^(pH - 6.8)));
          Real ZGSH=(-1) - 1*(10^(pH - 8.54)/(1 + 10^(pH - 8.54))) - 1*(10^(pH
               - 9.42)/(1 + 10^(pH - 9.42)));
          Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87))
            "this was not in the original implementatio by  machek!";
          //   Real Zim;//=-9.2 "Charge of ALL impermeable solutes";
          //   parameter Real Zim=-4.85 "Adjusted Charge of ALL impermeable solutes";
          parameter Real Zim=0;
          //-8.0452 "Charge of ALL impermeable solutes acc to MODEL";

          //osmolality
          Real fiHb=1 + 0.115*water_c[cont.Hb] + 0.0256*water_c[cont.Hb]^2
            "Osmotic coefficient of Hemoglobin (eq from Raftos et al)";
          Real fi[cont]={0.93,0.93,0.93,fiHb,1,1,1,1,1} "Osmotic coefficients";
          Real Osm=sum(water_c .* fi) + 0.93*(HCO3 + CO3) + Lactate +
              permeableParticles
            " = 284.084 osmotic coefficient from Wolf model V3.49";
          Real OsmTest[:]=water_c .* fi;
          //   //electric charge
          Real Z[cont]={1,1,-1,ZHb,ZDPG,ZATP,ZGSH,0,ZPi}
            "equivalent charges per mole";
          Real chargeTest[:]=water_c .* Z;
          Real EryConsIonCharge=sum(chargeTest[{cont.Hb,cont.DPG,cont.ATP,cont.GSH,
              cont.Pi}]);
          Real charge=(sum(water_c .* Z) - HCO3 - 2*CO3 - Lactate + Zim*Vew0/
              Vew)*Vew;

          // parameter Fraction Hct0 = 0.44;
          // parameter Volume totalBlodVolume;
          Real Vew0;
          constant Real wf0(unit="1") = 0.73
            "Fraction of water in the erythrocyte";

          Concentration HCO3=0.026*pCO2mmHg*10^(pH - 6.103)/few
            "HCO3 should be 16,8 acc to the model";
          Concentration CO3=HCO3*10^(pH - 10.2) "CO3 is 0,016";
          Real pH(start=7.22) = -log10(H);

          Real Vew(
            unit="l",
            min=0,
            max=10,
            start=1.43);
          Real pCO2mmHg(unit="1");
          Real H(
            start=10^(-7.2),
            min=0,
            max=1);
          Concentration Lactate "LACew =LACpw/rCl";
          parameter Concentration permeableParticles=10.64
            "glucose and urea concentration in PLasma water";

          Concentration HCO3ClDiff=HCO3/water_c[cont.Cl]*few;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end Erythrocyte;

        record Plasma
          import Modelica.SIunits.*;

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

          Concentration Ca=2.4;
          parameter Concentration Mg=0.8;
          //parameter Concentration Cl0=105;
          //  Real Cl_mass0(unit="mol") = Cl0/wf0*water_volume0;
          //  Real Cl_mass(unit="mol") = Cl/wf0*water_volume;
          parameter Concentration Pi=1.2;
          Concentration Alb=AlbPwGPerL/66.5 "mmol/Lpw";
          Real AlbPwGPerL=43*Vp0ByVp "g/lpw";
          // parameter Concentration im = 0; // original
          parameter Concentration im=11.87;
          // adjusted
          parameter Concentration Lac=1.5;
          Concentration SO4pw=0.33/fpw "mmol/Lpw";
          Concentration volume_c[cont]={Na,K,Ca,Mg,Cl,Pi,Alb/Vp0ByVp,im,Lac}
            "concentration in one liter";
          Concentration water_c[cont]=volume_c/fpw*Vp0ByVp
            "Actual concentration recalculated to water fraction in one liter (initially 0.96)";
          //   Real masses0[cont](each unit="mol") = volume_c ./ wf0 .* water_volume0
          //     "total mass of contents";
          Real CaPPBound=(3.98)/(3.98 + Hh)*water_c[cont.Ca] "1.257";
          Real MgPPBound=CaPPBound/2;
          Real ZCaBindPerAlb=CaPPBound/water_c[cont.Alb];
          Real ZClBindPerAlb=6.46/(6.46 + Hh)*6 + 4 "Anstey";
          Real CaIon=water_c[cont.Ca]*2 - CaPPBound "mEq/Lpw";
          Real MgIon=water_c[cont.Mg]*2 - MgPPBound "mEq/Lpw";
          //charge on inpermeable solutes
          Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87));
          Real ZFigge=(-10.65) - 16*(10^(pH - 7.418)/(1 + 10^(pH - 7.418)));
          Real ZAlbBnd=-ZClBindPerAlb + ZCaBindPerAlb + ZCaBindPerAlb/2;
          Real ZAlb=ZFigge + ZAlbBnd;
          Real Hh=H/fpw*1e8;

          //   parameter Real Zim ;//= -5.3 "Charge of ALL impermeable solutes";//test
          // Real Zim = -5.3 "Charge of ALL impermeable solutes";//test
          parameter Real Zim=-0.6349 "Charge of ALL impermeable solutes";

          Real fi[cont]={0.93,0.93,1,1,0.93,0.93,1,1,-9999};
          // Alb has osmotic coefficient of 1???
          Real CaOsm=water_c[cont.Ca] - CaPPBound*0.5;
          Real MgOsm=water_c[cont.Mg] - MgPPBound*0.5;
          Real OsmPart=SO4pw + CaOsm + MgOsm + sum(water_c[{cont.Na,cont.K,cont.Cl,
              cont.Lac}]) + water_c[cont.Pi] + HCO3 + CO3;
          Real Osm=OsmPart*0.93 + permeableParticles " = 284.084";

          Real Z[cont]={1,1,2,2,-1,ZPi,ZAlb,0,1};
          Real Chrgs[cont]=water_c .* Z;

          Real SID=sum(water_c[{cont.Na,cont.K}]) - water_c[cont.Cl] + 2*sum(
              water_c[{cont.Ca,cont.Mg}]) - CaPPBound - MgPPBound - water_c[
              cont.Lac];
          Real charge=Zim + SID + water_c[cont.Alb]*ZAlb + water_c[cont.Pi]*ZPi
               - HCO3 - 2*CO3 - SO4pw*2;
          //*water_volume;

          parameter Volume water_volume0=(1 - 0.44)*0.94*1;

          constant VolumeFraction wf0=0.94 "Fraction of water in the plasma";
          Concentration HCO3=0.0306*pCO2mmHg*10^(pH - 6.103)/fpw;
          Concentration CO3=HCO3*10^(pH - 10.2);
          Real pH(start=7.37) = -log10(H);

          Real Vp0ByVp;
          Real fpw "fraction plasma water - Vpw/Vp";
          Real pCO2mmHg(unit="1");
          Real H(
            start=10^(-7.37),
            min=0,
            max=1);
          Real Hw=H/fpw;
          Concentration permeableParticles=(5 + 5)/fpw
            "glucose and urea concentration in PLasma water";
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end Plasma;

        record Isf
          "ISF is considered as water compartment, that is each substance concentration is meant as in water"
          import Modelica.SIunits.*;

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
          Concentration MgIon "mEq/l";
          //Modelica.SIunits.Concentration CaPPBound "mEq/l";
          Concentration Cl;
          Real Vis;
          Real Vis0;
          parameter Real Alb(unit="g/l") = 43 "plasma alb conc g/l";
          Concentration SO4pw "mmol/Pw";

          Concentration albis=(Alb*0.33*Vis0)/(Vis - Vis0*0.25) "g/iw";
          Concentration albiw=albis/66.5
            "mmol/liw Alb conc is fixed as 1/3 of 43, no matter the conc in plasma";

          Real CaPPBound=(3.98)/(3.98 + Hh)*plasma_water_c[cont.Ca]
            "= 1.290 Meq/Lis";
          Real CaPPBoundErr=0
            "Error in VisSim diagram, should be CaPPBound instead";
          Real ZCaBindPerAlb=CaPPBound/albiw;
          Real ZClBindPerAlb=6.46/(6.46 + Hh)*6 + 4 "Anstey";

          Real transf[cont]={rClpwis,rClpwis,0,0,0,(1/rClpwis)^(-ZPi),0,0,1/
              rClpwis} "Only Na, K, Pi and Lac";

          Real SO4=(SO4pw*2)/(rClpwis^2) "concentration in one liter mEq/lw";
          // Concentration water_c[cont]=transf .* plasma_water_c;
          Concentration water_c[cont]=transf .* plasma_water_c
            "Actual concentration recalculated from plasma";

          Real ZPi=(-1) - 10^(pH - 6.87)/(1 + 10^(pH - 6.87));
          Real ZFigge=(-10.65) - 16*(10^(pH - 7.418)/(1 + 10^(pH - 7.418)));
          Real ZAlbBnd=-ZClBindPerAlb + ZCaBindPerAlb + ZCaBindPerAlb/2;
          Real ZAlb=ZFigge + ZAlbBnd;

          parameter Real Zim=3.057 "Charge of ALL impermeable solutes";

          Real OsmStep1=sum(water_c[{cont.Na,cont.K}]) + Cl;
          Real CaMgIonOsm=(CaIon - CaPPBoundErr + MgIon - CaPPBoundErr/2)/2;
          Real OsmStep2=OsmStep1 + CaMgIonOsm + sum(water_c[{cont.Pi,cont.Lac}])
               + HCO3 + CO3 + SO4/2;
          Real Osm=OsmStep2*0.93 + permeableParticles;

          Real SID=sum(water_c[{cont.Na,cont.K}]) - Cl + CaIon - CaPPBoundErr
               + (MgIon - 0.5*CaPPBoundErr) - water_c[cont.Lac];
          Real charge=Zim + SID + albiw*ZAlb + water_c[cont.Pi]*ZPi - HCO3 - 2*
              CO3 - SO4;
          //*water_volume;

          Concentration HCO3=0.0326*pCO2mmHg*10^(pH - 6.103);
          Concentration CO3=HCO3*10^(pH - 10.2);
          Real pH(start=7.408) = -log10(H);

          Real pCO2mmHg(unit="1");
          Real H(
            start=10^(-7.408),
            min=0,
            max=1);
          Real Hh=H*1e8;

          parameter Concentration permeableParticles=10.64
            "glucose and urea concentration in PLasma water";

          // Transcapillary pressures
          // wattenpaugh J trauma Inuj 1998
          Real AlbPwGPerL;
          Real ca=AlbPwGPerL*0.1 "Albumin concentration";
          Real cp=AlbPwGPerL*0.1*1.637 "protein concentration";
          Real COPPpl=PIEalb + PIEglob;
          Real PIEalb=ca/cp*(cp*2.8 + 0.18*cp^2 + 0.012*cp^3);
          Real PIEglob=(1 - ca/cp)*(0.9*cp + 0.12*cp^2 + 0.004*cp^3);
          Real cai=albis*0.1;
          Real COPalbis=cai*2.8 + 0.18*cai^2 + 0.012*cai^3;
          Real BloodyHell=(Vis - Vis0)/Vis0;
          Real Bh1=if BloodyHell < 0 then BloodyHell*125 else BloodyHell*243
            "Simplified bloody hell";
          //   Real pd = (28.48 - 0.99*(COPPpl - COPalbis))/19.3;
          Real pd=(38.36 - 0.99*(COPPpl - COPalbis) + PVb - Bh1)/19.3;
          Real OsmDiff=pd*2;
          Real Vbr=(Vb - Vb0)/Vb0 "VbRatio";
          Real PVb=if Vbr < -5 then -5*50 elseif Vbr > -5 and Vbr < 0 then Vbr*
              50 elseif Vbr > 0 and Vbr < 100 then Vbr*46.5 else 100*46.5
            "In original -5..0 and 0..100";
          Real Vb;
          Real Vb0;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end Isf;

        record Volumes
          type PatientCat = enumeration(
              adultM,
              adultF,
              childM,
              childF);
          parameter PatientCat pc=PatientCat.adultM;
          parameter Real Fcell=if pc == PatientCat.adultM or pc == PatientCat.adultF
               then 0.88 elseif pc == PatientCat.childM then 0.83 else 0.85;

          // the rest is only adult male
          parameter Real TBWm=2.447 - 0.09156*age + 0.3362*weight + 0.1074*
              height;
          //Real TBWf
          parameter Real LBM=0.407*weight + 0.267*height - 19.2;
          // Real LBMm, LBMf;
          parameter Real age=21;
          parameter Real weight=70;
          parameter Real height=175;

          parameter Real Ve0=0.048*LBM - 0.72;
          parameter Real Vp0=0.051*LBM + 0.27;
          parameter Real Vb0=(Ve0 + Vp0);
          parameter Real fH0=Ve0/(Ve0 + Vp0);
          parameter Real fHa0=fH0/Fcell;
          parameter Real Vis0=0.331*LBM - Vpw0;
          parameter Real TBW=2.447 - 0.09156*age + 0.3362*weight + 0.1074*
              height "Woodrow 2003 for male, there also exists female variant";
          parameter Real Vc0=TBW - Vis0 - Vpw0 - Vew0;

          parameter Real Vpw0=Vp0 - Vps;
          parameter Real Vps=0.0603*Vp0;
          parameter Real Vew0=Ve0*0.7317;
          parameter Real Ves=Ve0 - Vew0;
          Real Ve=Ves + Vew;
          Real Vb=Vpl + Ve;
          Real fH=Ve/Vb;

          parameter Real Vw0=Ve0 + Vp0 + Vis0 + Vc0;
          parameter Real V=Vw0 + Vadd;
          parameter Real Vadd=0;
          // Volumes
          Real Vpl=Vpw + Vps;
          Real fpw=Vpw/Vpl;
          Real few=Vew/Ve;
          Real Vpw=Vadd + Vpw0 + Vis0 + Vew0 + Vc0 - Vis - Vew - Vc " = 2.93";
          Real Vp0ByVp=Vp0/Vpl;
          Real Vew(start=Vew0,min=1e-3) "= 1.43";
          Real Vis(start=Vis0,min=1e-3) "= 15.59";
          Real Vc(start=Vc0,min=1e-3) "= 22.9";
          // Real Vew = 1.43;
          // Real Vis = 15.59;
          // Real Vc = 22.9;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end Volumes;

        record StrongIonMasses

          Real Ve0;
          Real Vp0;
          Real Vis0;
          Real Vc0;
          // any ion addition MUST be added here as well;
          constant Real Napl0=140;
          constant Real Nais0=141.5;
          constant Real Nac0=12;
          Real MNa=Napl0*Vp0 + Nais0*Vis0 + Nac0*Vc0;
          Real MNac=Nac0*Vc0;
          Real MNa_IP;

          constant Real Kpl0=4.7;
          constant Real Kis0=4.75;
          constant Real Kc0=139;
          Real MK=Kpl0*Vp0 + Kis0*Vis0 + Kc0*Vc0;
          Real MKc=Kc0*Vc0;
          Real MK_IP;

          constant Real Cle0=53.5;
          constant Real Clpl0=104;
          constant Real Clis0=116.5;
          constant Real Clc0=4;
          Real MCl=Cle0*Ve0 + Clpl0*Vp0 + Clis0*Vis0 + Clc0*Vc0;
          Real MClc=Clc0*Vc0;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end StrongIonMasses;

        record Cell
          import Modelica.SIunits.*;

          Real Kis;
          Real Nais;
          Real Clis;
          Real Vc0;
          Real Vc;
          Real His;
          Real pCO2mmHg(unit="1");

          Real Clc(
            start=3.9,
            min=1,
            max=100);

          Real rClisce(start=29.6) = Clis/Clc;
          Real Kc=Kis*rClisce;
          Real Nac=Nais*rClisce*0.0029;

          Concentration HCO3=0.029*pCO2mmHg*10^(pH - 6.103);
          Real pH(start=7.408) = -log10(H);

          Real H(
            start=10^(-6.935),
            min=1e-9,
            max=1) = His*0.1*rClisce "Why His*0.1? Is it the water content?";

          Real charge=-HCO3 + Kc - Clc + Nac + Zim*Im;
          Real Zim=-1.1728*10^(pH - 5.5)/(1 + 10^(pH - 5.5));
          Real Im=Vc0/Vc*123.045;

          Real OsmPart=HCO3 + Kc + Nac + Clc;
          Real Osm=Im + permeableParticles + 0.93*OsmPart;

          Real MNac=Nac*Vc;
          Real MClc=Clc*Vc;
          Real MKc=Kc*Vc;

          parameter Concentration permeableParticles=5.313 + UreaError
            "glucose and urea concentration in PLasma water";
          constant Real UreaError=0
            "Concentration of Urea in plasma water. Should be 5.31, is zero";

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end Cell;
      end Auxiliary;

      package Tests
        model E
          FullBloodAcidBase.Wolf.v351.Auxiliary.Erythrocyte ery;
          FullBloodAcidBase.Wolf.v351.Auxiliary.Volumes vols;
          // total mass of Cl mobile ion
          Real pCO2mmHg=45.1;
          // output: osmolarity and Cl
          parameter Real Hct=0.5;
          parameter Real O2s=0.70;

          parameter Real Hpw=4.602669e-8;
          parameter Real Clpw=110.6733;
          parameter Real EryOsm=286.03625;
          Real Cle=ery.water_c[ery.cont.Cl] "= 71.031139 @v3.49";
          Real Vew=vols.Vew "= 1.44400";
        equation
          vols.Vew = 1.44026;
          vols.Vis = 15.603374;
          vols.Vc = 22.784381;

          vols.Ve0/vols.Ve = ery.Ve0ByVeFraction;
          ery.fH = vols.fH;
          ery.few = vols.few;
          ery.Lactate = 1.083;
          //LACpw/rCl;
          ery.Vew0 = vols.Vew0;
          ery.Vew = vols.Vew;
          Cle/Clpw = Hpw/ery.H "Clpla , Hpl at standard V3.49";

          ery.water_c[ery.cont.Cl] = 75.1214;
          //1.0311;
          // ery.H = 6.37518e-8;

          // ery.
          //   // same osmolarities
          //  ery.Osm = EryOsm;
          //  ery.charge = 0;
          //   ery.Cl = ery.Cl0;
          //ery.Cl = 71.77 * ery.wf0;

          // Donnan equilibrium
          // pco2
          ery.pCO2mmHg = pCO2mmHg;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end E;

        model P
          Real Clpw(
            start=112,
            min=10,
            max=200) = pla.water_c[pla.cont.Cl] " = 109.694599";
          Real Vis=vols.Vis "=15.75516";

          FullBloodAcidBase.Wolf.v351.Auxiliary.Plasma pla;
          FullBloodAcidBase.Wolf.v351.Auxiliary.Volumes vols;
          FullBloodAcidBase.Wolf.v351.Auxiliary.StrongIonMasses sim;
          // total mass of Cl mobile ion
          Real pCO2mmHg=40;
          // Real rClpw_is = pla.volume_c[pla.cont.Cl]/Clis;
          Real rClpwis=Clpw/Clis;

          constant Real Clis=115.461419;
          constant Real MNac=272.17930;
          constant Real MKc=3178.3177;
          constant Real MClc=88.33349;
          constant Real MCle=102.5695;

          parameter Real EryOsm=284.084;
          Real MCl_IP=sim.MCl - MClc - MCle;
          //Real test = sim.MNa_IP /(vols.Vis * rClpw_is + vols.Vpw);
        equation
          vols.Vew = 1.44400;
          //  vols.Vis = 15.75516;
          vols.Vc = 22.64963;

          sim.Ve0 = vols.Ve0;
          sim.Vp0 = vols.Vp0;
          sim.Vis0 = vols.Vis0;
          sim.Vc0 = vols.Vc0;
          sim.MNa_IP = sim.MNa - MNac;
          sim.MK_IP = sim.MK - MKc;

          pla.water_c[pla.cont.Na] = sim.MNa_IP/(vols.Vis*rClpwis + vols.Vpw);
          pla.water_c[pla.cont.K] = sim.MK_IP/(vols.Vis*rClpwis + vols.Vpw);
          Clis = (MCl_IP - vols.Vpw*pla.water_c[pla.cont.Cl])/vols.Vis;
          pla.water_c[pla.cont.Cl] = 109.694599;

          // same osmolarities
          //  pla.Osm = EryOsm;
          pla.charge = 0;

          pla.Vp0ByVp = vols.Vp0ByVp;

          // Donnan equilibrium
          //  pla.pH = 7.410676;
          pla.fpw = vols.fpw;
          // pco2
          pla.pCO2mmHg = pCO2mmHg;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end P;

        model I
          FullBloodAcidBase.Wolf.v351.Auxiliary.Isf isf(permeableParticles=
                10.62755);
          FullBloodAcidBase.Wolf.v351.Auxiliary.Volumes vols;
          Real rClpwis=Clpw/isf.Cl " = 0.9500";
          parameter Real Hpw=4.128159e-8;
          parameter Real Clpw=109.694599;
          parameter Real Alb=43;
          constant Real EryOsm=284.084;
          constant Real MCl_IP=2148.68379;

          // Real Clis( start = 115) = MCl_IP - vols.Vpw*Clpw/vols.Vis " = 115.4614";
          Real Clis(start=115) " = 115.4614";
          Real Vis=vols.Vis "= 15.75516";
          Real oo=isf.Osm + isf.pd*2;
        equation
          vols.Vew = 1.44400;
          //  vols.Vis = 15.75516;
          vols.Vc = 22.64963;

          isf.plasma_water_c = {147.3244,4.9890,2.4977,0,0,1.2488,0,0,1.56017};

          isf.rClpwis = rClpwis;

          isf.CaIon = 3.76940*rClpwis^2;
          isf.MgIon = 1.05212*rClpwis^2;
          isf.Vis = vols.Vis;
          isf.Vis0 = vols.Vis0;
          isf.SO4pw = 0.702420/2;
          // The input here is in mmol/lpw to match plaswma variable, not in mEq/lpw as in the vissim model
          isf.pCO2mmHg = 40;
          isf.AlbPwGPerL = Alb*vols.Vp0ByVp;
          isf.Vb = vols.Vb;
          isf.Vb0 = vols.Vb0;
          isf.Cl = Clis;

          //  isf.permeableParticles = 10.62755 "5+5 / fpw";

          //  isf.Cl = 115.4614;
          //  isf.pH = 7.40649;

          isf.H = Hpw*rClpwis;
          isf.charge = 0;
          isf.Osm + isf.OsmDiff = EryOsm;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end I;

        model C
          import Modelica.SIunits.*;
          FullBloodAcidBase.Wolf.v351.Auxiliary.Cell c;
          Real Clc(
            start=3.9,
            min=1,
            max=100) = c.Clc;
        equation
          c.Kis = 4.73983;
          //4.73106;
          c.Nais = 139.9662;
          //141.474;
          c.Clis = 115.461419;
          //116.7911;
          c.Vc0 = 22.87176;
          //  c.Vc = 22.64963;//22.897789;
          c.His = 3.92197e-8;
          //3.91020e-8;
          c.pCO2mmHg = 40;

          //  c.rClisce = 29.6055;
          //  c.Clc = 3.961488;
          c.charge = 0;
          c.Osm = 282.20;

          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end C;
      end Tests;

      model EP

        Real Clpw(
          start=112,
          min=10,
          max=200) = pla.water_c[pla.cont.Cl] " = 109.694599";
        Real Cle(
          start=71,
          min=1,
          max=200) = ery.water_c[ery.cont.Cl] "= 71.031139 @v3.49";
        Real Vew=vols.Vew "= 1.44400";
        Real Vis=vols.Vis "=15.75516";

        Real pHe=ery.pH;
        Real pHp=pla.pH;

        FullBloodAcidBase.Wolf.v351.Auxiliary.Plasma pla;
        FullBloodAcidBase.Wolf.v351.Auxiliary.StrongIonMasses sim;
        FullBloodAcidBase.Wolf.v351.Auxiliary.Erythrocyte ery;
        FullBloodAcidBase.Wolf.v351.Auxiliary.Volumes vols;
        // total mass of Cl mobile ion
        Real pCO2mmHg=45.1;
        // output: osmolarity and Cl
        parameter Real Hct=0.50;
        parameter Real O2s=0.70;

        Real rClpwis=Clpw/Clis;

        constant Real Clis=116.513794;
        //115.461419;
        constant Real MNac=274.051633;
        constant Real MKc=3141.2286;
        constant Real MClc=91.4811;
        Real MCle=ery.water_c[ery.cont.Cl]*vols.Vew "= 102.5695";
        Real MCl_IP=sim.MCl - MClc - MCle;

      equation
        //  vols.Vew = 1.44400;
        //  vols.Vis = 15.75516;
        vols.Vc = 22.784381;

        vols.Ve0/vols.Ve = ery.Ve0ByVeFraction;
        ery.fH = vols.fH;
        ery.few = vols.few;
        ery.Lactate = pla.water_c[pla.cont.Lac]*Cle/Clpw;
        //LACpw/rCl;
        ery.Vew0 = vols.Vew0;
        ery.Vew = vols.Vew;

        Cle/Clpw = pla.Hw/ery.H "Clpla , Hpl at standard V3.49";
        ery.Osm = pla.Osm;
        ery.charge = 0;
        ery.pCO2mmHg = pCO2mmHg;

        sim.Ve0 = vols.Ve0;
        sim.Vp0 = vols.Vp0;
        sim.Vis0 = vols.Vis0;
        sim.Vc0 = vols.Vc0;
        sim.MNa_IP = sim.MNa - MNac;
        sim.MK_IP = sim.MK - MKc;

        pla.water_c[pla.cont.Na] = sim.MNa_IP/(vols.Vis*rClpwis + vols.Vpw);
        pla.water_c[pla.cont.K] = sim.MK_IP/(vols.Vis*rClpwis + vols.Vpw);
        Clis = (MCl_IP - vols.Vpw*pla.water_c[pla.cont.Cl])/vols.Vis;
        pla.water_c[pla.cont.Cl] = 110.6733;

        pla.charge = 0;
        pla.Vp0ByVp = vols.Vp0ByVp;
        pla.fpw = vols.fpw;
        pla.pCO2mmHg = pCO2mmHg;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end EP;

      model EPI

        Real Clis(
          start=115,
          min=10,
          max=200) = isf.Cl " = 115.4614";
        Real Clpw(
          start=112,
          min=10,
          max=200) = pla.water_c[pla.cont.Cl] " = 109.694599";
        Real Cle(
          start=71,
          min=1,
          max=200) = ery.water_c[ery.cont.Cl] "= 71.031139 @v3.49";

        Real Vew=vols.Vew "= 1.44400";
        Real Vis=vols.Vis "=15.75516";

        Real pHe=ery.pH;
        Real pHp=pla.pH;
        Real pHi=isf.pH;

        FullBloodAcidBase.Wolf.v351.Auxiliary.Plasma pla;
        FullBloodAcidBase.Wolf.v351.Auxiliary.StrongIonMasses sim;
        FullBloodAcidBase.Wolf.v351.Auxiliary.Erythrocyte ery;
        FullBloodAcidBase.Wolf.v351.Auxiliary.Isf isf(permeableParticles=
              10.62755);
        FullBloodAcidBase.Wolf.v351.Auxiliary.Volumes vols;
        // total mass of Cl mobile ion
        Real pCO2mmHg=40;
        // output: osmolarity and Cl
        parameter Real Hct=0.44;
        parameter Real O2s=0.75;

        constant Real MNac=272.17930;
        constant Real MKc=3178.3177;
        constant Real MClc=88.33349;
        Real MCle=ery.water_c[ery.cont.Cl]*vols.Vew "= 102.5695";
        Real MCl_IP=sim.MCl - MClc - MCle;
        Real rH=isf.H/pla.Hw;

        Real rClpwis=Clpw/Clis " = 0.9500";
        parameter Real Alb=43;

        // Real Clis( start = 115) = MCl_IP - vols.Vpw*Clpw/vols.Vis " = 115.4614";
      equation
        //  vols.Vew = 1.44400;
        //  vols.Vis = 15.75516;
        vols.Vc = 22.64963;

        vols.Ve0/vols.Ve = ery.Ve0ByVeFraction;
        ery.fH = vols.fH;
        ery.few = vols.few;
        ery.Lactate = pla.water_c[pla.cont.Lac]*Cle/Clpw;
        //LACpw/rCl;
        ery.Vew0 = vols.Vew0;
        ery.Vew = vols.Vew;

        Cle/Clpw = pla.Hw/ery.H "Clpla , Hpl at standard V3.49";
        ery.Osm = pla.Osm;
        ery.charge = 0;
        ery.pCO2mmHg = pCO2mmHg;

        sim.Ve0 = vols.Ve0;
        sim.Vp0 = vols.Vp0;
        sim.Vis0 = vols.Vis0;
        sim.Vc0 = vols.Vc0;
        sim.MNa_IP = sim.MNa - MNac;
        sim.MK_IP = sim.MK - MKc;

        pla.water_c[pla.cont.Na] = sim.MNa_IP/(vols.Vis*rClpwis + vols.Vpw);
        pla.water_c[pla.cont.K] = sim.MK_IP/(vols.Vis*rClpwis + vols.Vpw);
        Clis = (MCl_IP - vols.Vpw*pla.water_c[pla.cont.Cl])/vols.Vis;
        //pla.water_c[pla.cont.Cl] = 109.694599;

        pla.charge = 0;
        pla.Vp0ByVp = vols.Vp0ByVp;
        pla.fpw = vols.fpw;
        pla.pCO2mmHg = pCO2mmHg;

        isf.plasma_water_c = {147.3244,4.9890,2.4977,0,0,1.2488,0,0,1.56017};
        isf.rClpwis = rClpwis;
        isf.CaIon = 3.76940*rClpwis^2;
        isf.MgIon = 1.05212*rClpwis^2;
        isf.Vis = vols.Vis;
        isf.Vis0 = vols.Vis0;
        isf.SO4pw = pla.SO4pw;
        // The input here is in mmol/lpw to match plaswma variable, not in mEq/lpw as in the vissim model
        isf.pCO2mmHg = 40;
        isf.AlbPwGPerL = Alb*vols.Vp0ByVp;
        isf.Vb = vols.Vb;
        isf.Vb0 = vols.Vb0;

        isf.H = pla.Hw*rClpwis;
        isf.charge = 0;
        isf.Osm + isf.OsmDiff = pla.Osm;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end EPI;

      model EPIC
        Real Clc(
          start=3.9,
          min=1,
          max=100) = c.Clc;
        Real Clis(
          start=115,
          min=10,
          max=200) = isf.Cl " = 115.4614";
        Real Clpw(
          start=112,
          min=10,
          max=200) = pla.water_c[pla.cont.Cl] " = 109.694599";
        Real Cle(
          start=71,
          min=1,
          max=200) = ery.water_c[ery.cont.Cl] "= 71.031139 @v3.49";

        Real Vew=vols.Vew "= 1.44400";
        Real Vis=vols.Vis "=15.75516";

        Real pHe=ery.pH;
        Real pHp=pla.pH;
        Real pHi=isf.pH;
        Real pHc=c.pH;

        FullBloodAcidBase.Wolf.v351.Auxiliary.Plasma pla;
        FullBloodAcidBase.Wolf.v351.Auxiliary.StrongIonMasses sim;
        FullBloodAcidBase.Wolf.v351.Auxiliary.Erythrocyte ery;
        FullBloodAcidBase.Wolf.v351.Auxiliary.Isf isf(permeableParticles=
              10.62755);
        FullBloodAcidBase.Wolf.v351.Auxiliary.Volumes vols;
        FullBloodAcidBase.Wolf.v351.Auxiliary.Cell c;
        // total mass of Cl mobile ion
        Real pCO2mmHg=40;
        // output: osmolarity and Cl
        parameter Real Hct=0.44;
        parameter Real O2s=0.75;

        //   constant Real MNac=272.17930;
        //   constant Real MKc=3178.3177;
        //   constant Real MClc=88.33349;
        //   Real MNac=272.17930;
        //   Real MKc=3178.3177;
        //   Real MClc=88.33349;
        Real MNac=c.Nac*vols.Vc;
        Real MKc=c.Kc*vols.Vc;
        Real MClc=c.Clc*vols.Vc;
        Real MCle=ery.water_c[ery.cont.Cl]*vols.Vew "= 102.5695";
        Real MCl_IP=sim.MCl - MClc - MCle;
        Real rH=isf.H/pla.Hw;

        Real rClpwis=Clpw/Clis " = 0.9500";
        parameter Real Alb=43;

        // Real Clis( start = 115) = MCl_IP - vols.Vpw*Clpw/vols.Vis " = 115.4614";
      equation
        //  vols.Vew = 1.44400;
        //  vols.Vis = 15.75516;
        //  vols.Vc = 22.64963;

        vols.Ve0/vols.Ve = ery.Ve0ByVeFraction;
        ery.fH = vols.fH;
        ery.few = vols.few;
        ery.Lactate = pla.water_c[pla.cont.Lac]*Cle/Clpw;
        //LACpw/rCl;
        ery.Vew0 = vols.Vew0;
        ery.Vew = vols.Vew;

        Cle/Clpw = pla.Hw/ery.H "Clpla , Hpl at standard V3.49";
        ery.Osm = pla.Osm;
        ery.charge = 0;
        ery.pCO2mmHg = pCO2mmHg;

        sim.Ve0 = vols.Ve0;
        sim.Vp0 = vols.Vp0;
        sim.Vis0 = vols.Vis0;
        sim.Vc0 = vols.Vc0;
        sim.MNa_IP = sim.MNa - MNac;
        sim.MK_IP = sim.MK - MKc;

        pla.water_c[pla.cont.Na] = sim.MNa_IP/(vols.Vis*rClpwis + vols.Vpw);
        pla.water_c[pla.cont.K] = sim.MK_IP/(vols.Vis*rClpwis + vols.Vpw);
        Clis = (MCl_IP - vols.Vpw*pla.water_c[pla.cont.Cl])/vols.Vis;
        //pla.water_c[pla.cont.Cl] = 109.694599;

        pla.charge = 0;
        pla.Vp0ByVp = vols.Vp0ByVp;
        pla.fpw = vols.fpw;
        pla.pCO2mmHg = pCO2mmHg;

        isf.plasma_water_c = {147.3244,4.9890,2.4977,0,0,1.2488,0,0,1.56017};
        isf.rClpwis = rClpwis;
        isf.CaIon = 3.76940*rClpwis^2;
        isf.MgIon = 1.05212*rClpwis^2;
        isf.Vis = vols.Vis;
        isf.Vis0 = vols.Vis0;
        isf.SO4pw = pla.SO4pw;
        // The input here is in mmol/lpw to match plaswma variable, not in mEq/lpw as in the vissim model
        isf.pCO2mmHg = 40;
        isf.AlbPwGPerL = Alb*vols.Vp0ByVp;
        isf.Vb = vols.Vb;
        isf.Vb0 = vols.Vb0;

        isf.H = pla.Hw*rClpwis;
        isf.charge = 0;
        isf.Osm + isf.OsmDiff = pla.Osm;

        //    c.Kis = 4.73983;//4.73106;
        //    c.Nais = 139.9662;//141.474;
        //    c.Clis = 115.461419;//116.7911;
        //   c.Vc0 = 22.87176;
        //  c.Vc = 22.64963;//22.897789;
        //   c.His = 3.92197e-8;//3.91020e-8;
        //   c.pCO2mmHg = 40;

        //   c.Vc = vols.Vc;
        //  c.rClisce = 29.6055;
        //  c.Clc = 3.961488;
        //   c.charge = 0;
        //   c.Osm = 282.20;
        //

        c.Kis = isf.water_c[isf.cont.K];
        //4.73983;//4.73106;
        c.Nais = isf.water_c[isf.cont.Na];
        //139.9662;//141.474;
        c.Clis = isf.Cl;
        //115.461419;//116.7911;
        c.Vc0 = vols.Vc0;
        c.Vc = vols.Vc;
        c.His = isf.H;
        //3.92197e-8;//3.91020e-8;
        c.pCO2mmHg = 40;

        c.charge = 0;
        c.Osm = isf.Osm;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end EPIC;

      annotation (Documentation(info="<html>
<p>References:</p>
<p><br>Raftos, J. E., Bulliman, B. T., &amp; Kuchel, P. W. (1990). Evaluation of an electrochemical model of erythrocyte pH buffering using 31P nuclear magnetic resonance data. The Journal of General Physiology, 95(6), 1183&ndash;1204.</p>
<p>Wolf, M. B. (2013). Whole body acid-base and fluid-electrolyte balance: a mathematical model. American Journal of Physiology. Renal Physiology, 305(8), F1118&ndash;31.</p>
<p>Wolf, M. B. (2015). Comprehensive diagnosis of whole-body acid-base and fluid-electrolyte disorders using a mathematical model and whole-body base excess. Journal of Clinical Monitoring and Computing, 29(4), 475&ndash;490.</p>
<p>Wolf, M. B., &amp; DeLand, E. C. (2011). A comprehensive, computer-model-based approach for diagnosis and treatment of complex acid&ndash;base disorders in critically-ill patients. Journal of Clinical Monitoring and Computing, 25(6), 353&ndash;364.</p>
<p>Wolf, M. B., &amp; Deland, E. C. (2011). A mathematical model of blood-interstitial acid-base balance: application to dilution acidosis and acid-base status. Journal of Applied Physiology, 110(4), 988&ndash;1002.</p>
</html>"));
    end v351;
  end Wolf;

  annotation (uses(Physiolibrary(version="2.3.1"), Modelica(version="3.2.2")),
      Documentation(info="<html>
<p>Full-blood acid-base model package contains implementation of Figge-Fencl model, various formalizations of Siggaard-Andersen&apos;s nomogram and our proposed combination of both. In the Figures package, one could simulate all figures from the article &quot;<b>Modern and traditional acid-base approaches combined: a full blood model</b>&quot; , which this model is supplementing.</p>
<p>You can run all models with the green triangle. The rest are auxilliary subcomponents.</p>
<p>Physiolibrary version 2.3.1 is required to run the FullBloodSubmodelTest model. You can find it at physiolibrary.org</p>
<p>(C) Filip Ježek, Jiř&iacute; Kofr&aacute;nek under BSD 3-clause licence</p>
</html>"));
end FullBloodAcidBase;
