within WasteWater;
package BSM1 "Component models for the Activated Sludge Model No.1"
extends Modelica.Icons.Library;

model deni "ASM1 denitrification tank"
  //denitrification tank based on the ASM1 model

  extends WasteWater.Icons.deni;
  extends Interfaces.ASM1base;

  // tank specific parameters
  parameter Modelica.SIunits.Volume V=1000 "Volume of denitrification tank";

  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -10},{-90,10}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -10},{110,10}})));
  Interfaces.WWFlowAsm1out MeasurePort annotation (Placement(transformation(
          extent={{50,40},{60,50}})));
  Modelica.Blocks.Interfaces.RealInput T annotation (Placement(transformation(
          extent={{-110,30},{-90,50}})));
equation

  aeration = 0;
  // no aeration in this tank //

  // volume dependent dilution term of each concentration

  inputSi = (In.Si - Si)*In.Q/V;
  inputSs = (In.Ss - Ss)*In.Q/V;
  inputXi = (In.Xi - Xi)*In.Q/V;
  inputXs = (In.Xs - Xs)*In.Q/V;
  inputXbh = (In.Xbh - Xbh)*In.Q/V;
  inputXba = (In.Xba - Xba)*In.Q/V;
  inputXp = (In.Xp - Xp)*In.Q/V;
  inputSo = (In.So - So)*In.Q/V;
  inputSno = (In.Sno - Sno)*In.Q/V;
  inputSnh = (In.Snh - Snh)*In.Q/V;
  inputSnd = (In.Snd - Snd)*In.Q/V;
  inputXnd = (In.Xnd - Xnd)*In.Q/V;
  inputSalk = (In.Salk - Salk)*In.Q/V;

  annotation (
    Documentation(info="This component models the ASM1 processes and reactions taking place in an unaerated (denitrification) tank
of a wastewater treatment plant.

The InPort signal gives the tank temperature to the model.

Parameters:

    - all stoichiometric and kinetic parameters of the activated sludge model No.1 (ASM1)
  V - volume of the tank [m3]
"));
end deni;

model nitri "ASM1 nitrification tank"
  // nitrification (aerated) tank, based on the ASM1 model

  extends WasteWater.Icons.nitri;
  extends Interfaces.ASM1base;

  // tank specific parameters
  parameter Modelica.SIunits.Volume V=1333 "Volume of nitrification tank";

  // aeration system dependent parameters
  parameter Real alpha=0.7 "Oxygen transfer factor";
  parameter Modelica.SIunits.Length de=4.5 "depth of aeration";
  parameter Real R_air=23.5 "specific oxygen feed factor [gO2/(m^3*m)]";

  WWU.MassConcentration So_sat "Dissolved oxygen saturation";

  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -10},{-90,10}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -10},{110,10}})));
  Interfaces.WWFlowAsm1out MeasurePort annotation (Placement(transformation(
          extent={{50,40},{60,50}})));
  Modelica.Blocks.Interfaces.RealInput T annotation (Placement(transformation(
          extent={{-110,30},{-90,50}})));
  Interfaces.AirFlow AirIn annotation (Placement(transformation(extent={{-5,
            -103},{5,-93}})));
  parameter Real Kla = 240;
equation

  // Temperature dependent oxygen saturation by Simba
  // So_sat =13.89 + (-0.3825 + (0.007311 - 0.00006588*T)*T)*T;
  So_sat = 8;

  // extends the Oxygen differential equation by an aeration term
  // aeration [mgO2/l]; AirIn.Q_air needs to be in
  // Simulationtimeunit [m3*day^-1]
  // aeration = (alpha*(So_sat - So)/So_sat*AirIn.Q_air*R_air*de)/V;
  aeration = Kla * (So_sat - So);

  // volume dependent dilution term of each concentration

  inputSi = (In.Si - Si)*In.Q/V;
  inputSs = (In.Ss - Ss)*In.Q/V;
  inputXi = (In.Xi - Xi)*In.Q/V;
  inputXs = (In.Xs - Xs)*In.Q/V;
  inputXbh = (In.Xbh - Xbh)*In.Q/V;
  inputXba = (In.Xba - Xba)*In.Q/V;
  inputXp = (In.Xp - Xp)*In.Q/V;
  inputSo = (In.So - So)*In.Q/V;
  inputSno = (In.Sno - Sno)*In.Q/V;
  inputSnh = (In.Snh - Snh)*In.Q/V;
  inputSnd = (In.Snd - Snd)*In.Q/V;
  inputXnd = (In.Xnd - Xnd)*In.Q/V;
  inputSalk = (In.Salk - Salk)*In.Q/V;

  annotation (
    Documentation(info="This component models the ASM1 processes and reactions taking place in an aerated (nitrification) tank
of a wastewater treatment plant.

The InPort signal gives the tank temperature to the model.

Parameters:

        - all soichiometric and kinetic parameters of the activated sludge model No.1 (ASM1)
  V     - volume of the tank [m3]
  alpha - oxygen transfer factor
  de    - depth of the aeration system [m]
  R_air - specific oxygen feed factor [g O2/(m3*m)]
"), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
        graphics));
end nitri;

model SecClarModTakacs "Secondary Clarifier ASM1 Model based on Takacs"

  extends WasteWater.Icons.SecClar;
  extends SecClar.Takacs.Interfaces.ratios;
  package SCP = SecClar.Takacs;
  import SI = Modelica.SIunits;
  package WI = WasteWater.BSM1.Interfaces;
  package WWU = WasteWater.WasteWaterUnits;
  parameter SI.Length hsc=4.0 "height of secondary clarifier";
  parameter Integer n=10 "number of layers of SC model";
  parameter SI.Length zm=hsc/(1.0*n) "height of m-th secondary clarifier layer";
  parameter SI.Area Asc=1500.0 "area of secondary clarifier";
  parameter WWU.MassConcentration Xt=3000.0 "threshold for X";

  // total sludge concentration in clarifier feed
  WWU.MassConcentration Xf;

  WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},{
            -90,24}})));
  WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,47},
            {112,67}})));
  WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
            -106},{-20,-86}})));
  WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,-106},
            {40,-86}})));

  // layers 1 to 10
  SCP.bottom_layer S1(
    zm=zm,
    Asc=Asc,
    Xf=Xf,
    rXs=rXs,
    rXbh=rXbh,
    rXba=rXba,
    rXp=rXp,
    rXi=rXi,
    rXnd=rXnd) annotation (Placement(transformation(extent={{-35,-93},{35,-78}})));
  SCP.lower_layer S2(
    zm=zm,
    Asc=Asc,
    Xf=Xf) annotation (Placement(transformation(extent={{-35,-74},{35,-59}})));
  SCP.lower_layer S3(
    zm=zm,
    Asc=Asc,
    Xf=Xf) annotation (Placement(transformation(extent={{-35,-55},{35,-40}})));
  SCP.lower_layer S4(
    zm=zm,
    Asc=Asc,
    Xf=Xf) annotation (Placement(transformation(extent={{-35,-36},{35,-21}})));
  SCP.lower_layer S5(
    zm=zm,
    Asc=Asc,
    Xf=Xf) annotation (Placement(transformation(extent={{-35,-17},{35,-2}})));
  SCP.feed_layer S6(
    zm=zm,
    Asc=Asc,
    Xf=Xf) annotation (Placement(transformation(extent={{-35,2},{35,17}})));
  SCP.upper_layer S7(
    zm=zm,
    Asc=Asc,
    Xf=Xf,
    Xt=Xt) annotation (Placement(transformation(extent={{-35,21},{35,36}})));
  SCP.upper_layer S8(
    zm=zm,
    Asc=Asc,
    Xt=Xt,
    Xf=Xf) annotation (Placement(transformation(extent={{-35,40},{35,55}})));
  SCP.upper_layer S9(
    zm=zm,
    Asc=Asc,
    Xf=Xf,
    Xt=Xt) annotation (Placement(transformation(extent={{-35,59},{35,74}})));
  SCP.top_layer S10(
    zm=zm,
    Asc=Asc,
    Xf=Xf,
    Xt=Xt,
    rXs=rXs,
    rXbh=rXbh,
    rXba=rXba,
    rXp=rXp,
    rXi=rXi,
    rXnd=rXnd) annotation (Placement(transformation(extent={{-35,78},{35,93}})));
equation

  connect(S1.Up, S2.Dn) annotation (Line(points={{-2.22045e-15,-78},{
          -2.22045e-15,-74}}));
  connect(S2.Up, S3.Dn) annotation (Line(points={{-2.22045e-15,-59},{
          -2.22045e-15,-55}}));
  connect(S3.Up, S4.Dn) annotation (Line(points={{-2.22045e-15,-40},{
          -2.22045e-15,-36}}));
  connect(S5.Up, S6.Dn) annotation (Line(points={{-2.22045e-15,-2},{
          -2.22045e-15,2}}));
  connect(S6.Up, S7.Dn) annotation (Line(points={{-2.22045e-15,17},{
          -2.22045e-15,21}}));
  connect(S7.Up, S8.Dn) annotation (Line(points={{-2.22045e-15,36},{
          -2.22045e-15,40}}));
  connect(S9.Up, S10.Dn) annotation (Line(points={{-2.22045e-15,74},{
          -2.22045e-15,78}}));
  connect(S4.Up, S5.Dn) annotation (Line(points={{-2.22045e-15,-21},{
          -2.22045e-15,-17}}));
  connect(S8.Up, S9.Dn) annotation (Line(points={{-2.22045e-15,55},{
          -2.22045e-15,59}}));
  connect(Feed, S6.In) annotation (Line(points={{-100,10},{-67.5,10},{-67.5,9.8},
          {-35,9.8}}));
  connect(S1.PQw, Waste) annotation (Line(points={{17.5,-93},{17.5,-100},{30,
          -100}}));
  connect(S10.Out, Effluent) annotation (Line(points={{35,85.5},{67.5,85.5},{
          67.5,57},{100,57}}));
  connect(S1.PQr, Return) annotation (Line(points={{-21,-93},{-21,-100},{-30,
          -100}}));

  // total sludge concentration in clarifier feed
  Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

  // ratios of solid components
  rXs = Feed.Xs/Xf;
  rXbh = Feed.Xbh/Xf;
  rXba = Feed.Xba/Xf;
  rXp = Feed.Xp/Xf;
  rXi = Feed.Xi/Xf;
  rXnd = Feed.Xnd/Xf;

  annotation (
    Documentation(info="This component models an ASM1 10 - layer secondary clarifier model with 4 layers above the feed_layer (including top_layer)
and 5 layers below the feed_layer (including bottom_layer) based on Takac`s theory.

Parameters:
  hsc -  height of clarifier [m]
  n   -  number of layers
  Asc -  surface area of sec. clar. [m2]
  Xt  -  threshold value for Xtss [mg/l]

"));
end SecClarModTakacs;

model blower "Blower for the aeration of the nitrification tanks"

  extends WasteWater.Icons.blower;
  package WWU = WasteWater.WasteWaterUnits;

  parameter WWU.VolumeFlowRate Q_max=20000 "maximum blower capacity";
  parameter WWU.VolumeFlowRate Q_min=0.0 "minimum blower capacity";
  Real H;
  // this is just a help variable to reduce expressions

  Interfaces.AirFlow AirOut annotation (Placement(transformation(extent={{-20,
            90},{0,110}})));
  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(
        origin={98,-30},
        extent={{-10,-10},{10,10}},
        rotation=180)));
equation

  H =0.5*(-Q_min + Q_max) + u*0.5*(-Q_min + Q_max) + Q_min;
  AirOut.Q_air = -(if H > Q_max then Q_max else if H < Q_min then Q_min else H);

  annotation (
    Documentation(info="This component models a blower of a wastewater treatment plant which generates an airflow that is needed
for the nitrification.
The blower is connected to the nitrification tank.
The airflow is controlled by a signal u (-1 <= u <= 1).

Parameter:

  Qmax - maximum blower capacity [m3 Air/d], this is produced when the control signal u is 1 or greater.
  Qmin - minimum blower capacity [m3 Air/d], this is produced when the control signal u is -1 or below.

"));
end blower;

model pump "ASM1 wastewater pump"

  extends WasteWater.Icons.pump;
  package WWU = WasteWater.WasteWaterUnits;

  parameter WWU.VolumeFlowRate Q_min=0.0 "minimum pump capacity";
  parameter WWU.VolumeFlowRate Q_max=20000 "maximum pump capacity";
  Real H;
  // this is just a help variable to reduce expressions

  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -43},{-90,-23}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            18},{110,38}})));
  Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
          extent={{-99,15},{-79,35}})));
equation

  H =0.5*(-Q_min + Q_max) + u*0.5*(-Q_min + Q_max) + Q_min;
  Out.Q = -(if H > Q_max then Q_max else if H < Q_min then Q_min else H);

  Out.Q + In.Q = 0;
  Out.Si = In.Si;
  Out.Ss = In.Ss;
  Out.Xi = In.Xi;
  Out.Xs = In.Xs;
  Out.Xbh = In.Xbh;
  Out.Xba = In.Xba;
  Out.Xp = In.Xp;
  Out.So = In.So;
  Out.Sno = In.Sno;
  Out.Snh = In.Snh;
  Out.Snd = In.Snd;
  Out.Xnd = In.Xnd;
  Out.Salk = In.Salk;

  annotation (
    Documentation(info="This component models an ASM1 wastewater pump. It generates a wastewater flow
that is controlled by the signal u (-1 <= u <=1).

Parameter:

  Qmax - maximum pump capacity [m3/d], this is produced when the control signal u is 1 or greater.
  Qmin - minimum pump capacity [m3/d], this is produced when the control signal u is -1 or below.

"));
end pump;

model FlowSource "Flowsource"

  extends WasteWater.Icons.FlowSource;
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{88,
            -80},{108,-60}})));
  Modelica.Blocks.Interfaces.RealInput data
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
equation

  Out.Q =-data;

  annotation (
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Ellipse(
          extent={{-54,54},{56,-54}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-8,-54},{-14,-52},{-24,-48},{-32,-44},{-36,-40},{-42,-34},{
              -48,-26},{-50,-20},{52,-20},{50,-26},{46,-32},{42,-36},{38,-40},{
              34,-44},{30,-46},{26,-48},{22,-50},{16,-52},{10,-54},{4,-54},{0,
              -54},{-8,-54}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,95,191},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-54,54},{56,-54}},
          lineColor={0,0,0},
          lineThickness=0.5),
        Rectangle(
          extent={{-4,-52},{4,-74}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,95,191},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{4,-74},{88,-68}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,95,191},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-4,-54},{-4,-74},{88,-74}},
          thickness=0.5),
        Line(
          points={{4,-54},{4,-68},{88,-68}},
          thickness=0.5)}),
    Documentation(info="This component is used to feed an ASM1 wwtp model with flow data from measurement
when e.g. concentration is measured after the primary clarifier.

The dimension of InPort is 1.

  1 volumeflowrate Q of incoming wastewater [m3/d]"));
end FlowSource;

model WWSource "Wastewater source"

  extends WasteWater.Icons.WWSource;
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{88,
            -80},{108,-60}})));
  Modelica.Blocks.Interfaces.RealInput data[14]
    annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
equation

  Out.Q =-18446;
  Out.Si =30.0;
  Out.Ss =69.50;
  Out.Xi =51.20;
  Out.Xs =202.32;
  Out.Xbh =28.17;
  Out.Xba =0;
  Out.Xp =0;
  Out.So =0;
  Out.Sno =0;
  Out.Snh =31.56;
  Out.Snd =6.95;
  Out.Xnd =10.59;
  Out.Salk =7.0;

  annotation (
    Documentation(info="This component provides all ASM1 data at the influent of a wastewater treatment plant.
The dimension of InPort is 14.

  1 volumeflowrate Q of incoming wastewater [m3/d]
  2 Si  [g COD/m3]
  3 Ss  [g COD/m3]
  4 Xi  [g COD/m3]
  5 Xs  [g COD/m3]
  6 Xbh [g COD/m3]
  7 Xba [g COD/m3]
  8 Xp  [g COD/m3]
  9 So  [g O2/m3]
 10 Sno [g N/m3]
 11 Snh [g N/m3]
 12 Snd [g N/m3]
 13 Xnd [g N/m3]
 14 Salk[mmol/l]
"));
end WWSource;

model EffluentSink "Receiving water (river)"
  // only for graphical termination in diagram layer, no equation needed

  extends WasteWater.Icons.EffluentSink;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            10},{-90,30}})));
  annotation (
    Documentation(info="This component terminates an ASM1 wastewater treatment plant model e.g. the wastewater flow to the receiving water.
"));
end EffluentSink;

model SludgeSink "Wastesludge sink"
  // only for graphical termination in diagram layer, no equation needed

  extends WasteWater.Icons.SludgeSink;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -22},{-90,-2}})));
  annotation (
    Documentation(info="This component terminates the waste sludge stream of an ASM1 wastewater treatment plant model.
Storage or further sludge treatment is not jet considered."));
end SludgeSink;

model ControlledDivider2 "Controlled flow divider"
  // divides one flow of wastewater into 2 Flows controlled by the

  // input signal u; u=1 means Out1.Q=In.Q and u=0 means Out2.Q=In.Q

  extends WasteWater.Icons.ControlledDivider2;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-111,
            -7},{-91,13}})));
  Interfaces.WWFlowAsm1out Out1 annotation (Placement(transformation(extent={{
            90,15},{110,35}})));
  Interfaces.WWFlowAsm1out Out2 annotation (Placement(transformation(extent={{
            90,-25},{110,-5}})));
  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(
        origin={0,-60},
        extent={{-10,-10},{10,10}},
        rotation=90)));
equation

  Out1.Q =-In.Q*u;
  Out2.Q =-In.Q*(1 - u);

  Out1.Si = In.Si;
  Out1.Ss = In.Ss;
  Out1.Xi = In.Xi;
  Out1.Xs = In.Xs;
  Out1.Xbh = In.Xbh;
  Out1.Xba = In.Xba;
  Out1.Xp = In.Xp;
  Out1.So = In.So;
  Out1.Sno = In.Sno;
  Out1.Snh = In.Snh;
  Out1.Snd = In.Snd;
  Out1.Xnd = In.Xnd;
  Out1.Salk = In.Salk;

  Out2.Si = In.Si;
  Out2.Ss = In.Ss;
  Out2.Xi = In.Xi;
  Out2.Xs = In.Xs;
  Out2.Xbh = In.Xbh;
  Out2.Xba = In.Xba;
  Out2.Xp = In.Xp;
  Out2.So = In.So;
  Out2.Sno = In.Sno;
  Out2.Snh = In.Snh;
  Out2.Snd = In.Snd;
  Out2.Xnd = In.Xnd;
  Out2.Salk = In.Salk;

  annotation (
    Documentation(info="This component divides one wastewater flow (ASM1) into two flows which are controlled by the signal u (0...1).
Is u.signal=1, the flow goes to output 1 (Out1) and is u.signal=0, the flow goes to output 2 (Out2).
The concentrations of the outport-flows are equal to the concentration at inport."));
end ControlledDivider2;

model divider2 "Flowdivider"

    // divides one flow of wastewater into 2 Flows; one amount needs to be specified

  extends WasteWater.Icons.divider2;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,-9},
              {-90,11}})));
  Interfaces.WWFlowAsm1out Out1 annotation (Placement(transformation(extent={{90,15},
              {110,35}})));
  Interfaces.WWFlowAsm1out Out2 annotation (Placement(transformation(extent={{90,-26},
              {110,-6}})));
equation

  In.Q + Out1.Q + Out2.Q = 0;

  Out1.Si = In.Si;
  Out1.Ss = In.Ss;
  Out1.Xi = In.Xi;
  Out1.Xs = In.Xs;
  Out1.Xbh = In.Xbh;
  Out1.Xba = In.Xba;
  Out1.Xp = In.Xp;
  Out1.So = In.So;
  Out1.Sno = In.Sno;
  Out1.Snh = In.Snh;
  Out1.Snd = In.Snd;
  Out1.Xnd = In.Xnd;
  Out1.Salk = In.Salk;

  Out2.Si = In.Si;
  Out2.Ss = In.Ss;
  Out2.Xi = In.Xi;
  Out2.Xs = In.Xs;
  Out2.Xbh = In.Xbh;
  Out2.Xba = In.Xba;
  Out2.Xp = In.Xp;
  Out2.So = In.So;
  Out2.Sno = In.Sno;
  Out2.Snh = In.Snh;
  Out2.Snd = In.Snd;
  Out2.Xnd = In.Xnd;
  Out2.Salk = In.Salk;

  annotation (
    Documentation(info=
          "This component divides one ASM1 wastewater flow into two ASM1 wastewater flows."), Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics));
end divider2;

model mixer2 "Mixer of two ASM1 characterised flows"

  extends WasteWater.Icons.mixer2;
  Interfaces.WWFlowAsm1in In1 annotation (Placement(transformation(extent={{
            -110,15},{-90,35}})));
  Interfaces.WWFlowAsm1in In2 annotation (Placement(transformation(extent={{
            -110,-25},{-90,-5}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -5},{110,15}})));
equation

  In1.Q + In2.Q + Out.Q = 0;
  Out.Si = (In1.Si*In1.Q + In2.Si*In2.Q)/(In1.Q + In2.Q);
  Out.Ss = (In1.Ss*In1.Q + In2.Ss*In2.Q)/(In1.Q + In2.Q);
  Out.Xi = (In1.Xi*In1.Q + In2.Xi*In2.Q)/(In1.Q + In2.Q);
  Out.Xs = (In1.Xs*In1.Q + In2.Xs*In2.Q)/(In1.Q + In2.Q);
  Out.Xbh = (In1.Xbh*In1.Q + In2.Xbh*In2.Q)/(In1.Q + In2.Q);
  Out.Xba = (In1.Xba*In1.Q + In2.Xba*In2.Q)/(In1.Q + In2.Q);
  Out.Xp = (In1.Xp*In1.Q + In2.Xp*In2.Q)/(In1.Q + In2.Q);
  Out.So = (In1.So*In1.Q + In2.So*In2.Q)/(In1.Q + In2.Q);
  Out.Sno = (In1.Sno*In1.Q + In2.Sno*In2.Q)/(In1.Q + In2.Q);
  Out.Snh = (In1.Snh*In1.Q + In2.Snh*In2.Q)/(In1.Q + In2.Q);
  Out.Snd = (In1.Snd*In1.Q + In2.Snd*In2.Q)/(In1.Q + In2.Q);
  Out.Xnd = (In1.Xnd*In1.Q + In2.Xnd*In2.Q)/(In1.Q + In2.Q);
  Out.Salk = (In1.Salk*In1.Q + In2.Salk*In2.Q)/(In1.Q + In2.Q);

  annotation (
    Documentation(info=
          "This component mixes two flows of wastewater (ASM1) of different concentration and different amount."));
end mixer2;

model mixer3 "Mixer of 3 ASM1 characterised flows"

  extends WasteWater.Icons.mixer3;
  Interfaces.WWFlowAsm1in In1 annotation (Placement(transformation(extent={{
            -110,25},{-90,45}})));
  Interfaces.WWFlowAsm1in In2 annotation (Placement(transformation(extent={{
            -110,-15},{-90,5}})));
  Interfaces.WWFlowAsm1in In3 annotation (Placement(transformation(extent={{
            -110,-55},{-90,-35}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -14},{110,6}})));
equation

  In1.Q + In2.Q + In3.Q + Out.Q = 0;
  Out.Si = (In1.Si*In1.Q + In2.Si*In2.Q + In3.Si*In3.Q)/(In1.Q + In2.Q + In3.Q);
  Out.Ss = (In1.Ss*In1.Q + In2.Ss*In2.Q + In3.Ss*In3.Q)/(In1.Q + In2.Q + In3.Q);
  Out.Xi = (In1.Xi*In1.Q + In2.Xi*In2.Q + In3.Xi*In3.Q)/(In1.Q + In2.Q + In3.Q);
  Out.Xs = (In1.Xs*In1.Q + In2.Xs*In2.Q + In3.Xs*In3.Q)/(In1.Q + In2.Q + In3.Q);
  Out.Xbh = (In1.Xbh*In1.Q + In2.Xbh*In2.Q + In3.Xbh*In3.Q)/(In1.Q + In2.Q +
    In3.Q);
  Out.Xba = (In1.Xba*In1.Q + In2.Xba*In2.Q + In3.Xba*In3.Q)/(In1.Q + In2.Q +
    In3.Q);
  Out.Xp = (In1.Xp*In1.Q + In2.Xp*In2.Q + In3.Xp*In3.Q)/(In1.Q + In2.Q + In3.Q);
  Out.So = (In1.So*In1.Q + In2.So*In2.Q + In3.So*In3.Q)/(In1.Q + In2.Q + In3.Q);
  Out.Sno = (In1.Sno*In1.Q + In2.Sno*In2.Q + In3.Sno*In3.Q)/(In1.Q + In2.Q +
    In3.Q);
  Out.Snh = (In1.Snh*In1.Q + In2.Snh*In2.Q + In3.Snh*In3.Q)/(In1.Q + In2.Q +
    In3.Q);
  Out.Snd = (In1.Snd*In1.Q + In2.Snd*In2.Q + In3.Snd*In3.Q)/(In1.Q + In2.Q +
    In3.Q);
  Out.Xnd = (In1.Xnd*In1.Q + In2.Xnd*In2.Q + In3.Xnd*In3.Q)/(In1.Q + In2.Q +
    In3.Q);
  Out.Salk = (In1.Salk*In1.Q + In2.Salk*In2.Q + In3.Salk*In3.Q)/(In1.Q + In2.Q
     + In3.Q);

  annotation (
    Documentation(info=
          "This component mixes 3 flows of wastewater (ASM1) of different concentration and different amount."));
end mixer3;

model sensor_COD "Ideal sensor to measure chemical oxygen demand (COD)"

  extends WasteWater.Icons.sensor_COD;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput COD annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));
equation

  In.Q = 0.0;
  COD = In.Si + In.Ss + In.Xi + In.Xs + In.Xbh + In.Xba + In.Xp;

  annotation (
    Documentation(info="This component measures the chemical oxygen demand (COD) concentration [g/m3]
of ASM1 wastewater and provides the result as output signal (to be
further processed with blocks of the Modelica.Blocks library).
"));
end sensor_COD;

model sensor_NH "Ideal sensor to measure ammonium nitrogen"

  extends WasteWater.Icons.sensor_NH;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput Snh annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));
equation

  In.Q = 0;
  Snh = In.Snh;

  annotation (
    Documentation(info="This component measures the ammonium nitrogen concentration [g/m3]
of ASM1 wastewater and provides the result as output signal (to be
further processed with blocks of the Modelica.Blocks library).
"));
end sensor_NH;

model sensor_NO "Ideal sensor to measure nitrate nitrogen"

  extends WasteWater.Icons.sensor_NO;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput Sno annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));
equation

  In.Q = 0;
  Sno = In.Sno;

  annotation (
    Documentation(info="This component measures the nitrate nitrogen concentration [g/m3]
of ASM1 wastewater and provides the result as output signal (to be
further processed with blocks of the Modelica.Blocks library).
"));
end sensor_NO;

model sensor_O2 "Ideal sensor to measure dissolved oxygen concentration"

  extends WasteWater.Icons.sensor_O2;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput So annotation (Placement(transformation(
          extent={{88,-10},{108,10}})));
equation

  In.Q = 0;
  So = In.So;

  annotation (
    Documentation(info="This component measures the dissolved oxygen concentration [g/m3]
of ASM1 wastewater and provides the result as output signal (to be
further processed with blocks of the Modelica.Blocks library).
"), Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Ellipse(
          extent={{-50,50},{50,-50}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={223,223,159},
          fillPattern=FillPattern.Solid),
        Line(
          points={{0,50},{0,38}},
          thickness=0.5),
        Line(
          points={{-50,0},{38,0}},
          thickness=0.5),
        Line(
          points={{50,0},{38,0}},
          thickness=0.5),
        Line(
          points={{-36,34},{-28,26}},
          thickness=0.5),
        Line(
          points={{34,36},{26,28}},
          thickness=0.5),
        Line(
          points={{0,0},{26,28}},
          thickness=0.5),
        Polygon(
          points={{30,32},{10,24},{24,12},{30,32}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Text(extent={{-36,-10},{36,-32}}, textString=
                                              "O2"),
        Line(
          points={{0,-50},{0,-90}},
          thickness=0.5),
        Line(points={{50,0},{88,0}}),
        Text(extent={{-80,100},{80,60}}, textString=
                                             "%name")}));
end sensor_O2;

model sensor_Q
    "Ideal sensor to measure the flow rate of an ASM1 wastewater stream"

  extends WasteWater.Icons.sensor_Q;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -10},{-90,10}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -10},{110,10}})));
  Modelica.Blocks.Interfaces.RealOutput Q
    annotation (Placement(transformation(
        origin={0,-98},
        extent={{-10,-10},{10,10}},
        rotation=270)));
equation

  In.Q + Out.Q = 0;
  Q = In.Q;
  // eventually abs(In.Q) to be shure to have pos. signal
  In.Si = Out.Si;
  In.Ss = Out.Ss;
  In.Xi = Out.Xi;
  In.Xs = Out.Xs;
  In.Xbh = Out.Xbh;
  In.Xba = Out.Xba;
  In.Xp = Out.Xp;
  In.So = Out.So;
  In.Sno = Out.Sno;
  In.Snh = Out.Snh;
  In.Snd = Out.Snd;
  In.Xnd = Out.Xnd;
  In.Salk = Out.Salk;

  annotation (
    Documentation(info="This component measures the flow of an ASM1 wastewater stream and provides
the result as output signal (to be further processed with blocks of
the Modelica.Blocks library).
"), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}),
        graphics));
end sensor_Q;

model sensor_TKN "Ideal TKN and total nitrogen sensor"

  extends WasteWater.Icons.sensor_TKN;
  extends Interfaces.stoichiometry;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput TKN[2]
    annotation (Placement(transformation(extent={{88,-10},{108,10}})));
equation

  In.Q = 0.0;
  TKN[1] = In.Snh + In.Snd + In.Xnd + i_xb*(In.Xbh + In.Xba) + i_xp*(In.Xp + In.Xi);
  TKN[2] = TKN[1] + In.Sno;

  annotation (
    Documentation(info="This component measures the Total Kjeldal Nitrogen (TKN) and the
total nitrogen (N_total) concentration [g/m3] of ASM1 wastewater
and provides the result as output signal (to be further processed
with blocks of the Modelica.Blocks library).

signal[1] - TKN
signal[2] - N_total
"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics));
end sensor_TKN;

model sensor_TSS
    "Ideal sensor to measure total suspended solids concentration (ASM1)"

  extends WasteWater.Icons.sensor_TSS;
  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput TSS annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));
equation

  In.Q = 0;

  TSS = 0.75*(In.Xs + In.Xbh + In.Xba + In.Xp + In.Xi);
  // the factor 0.75 needs to be adapted due to plant dependency
  // 0.75 is from the COST Benchmark configuration

  annotation (
    Documentation(info="This component measures the total suspended solids concentration [g/m3]
of ASM1 wastewater and provides the result as output signal (to be
further processed with blocks of the Modelica.Blocks library).
"));
end sensor_TSS;

  package Examples
    "Demonstration examples of the components of the ASM1 library"

    extends Modelica.Icons.Library;

    class SmallPlant "Small WWTP Configuration"
      import WasteWater;
      extends Modelica.Icons.Example;

      //Q_air=12100.99290780142 is equal to a Kla of 3.5 h^-1 from COST benchmark
      //Q_air=34574.2654508612 is equal to a Kla of 10 h^-1 from COST benchmark

      WasteWater.BSM1.EffluentSink Effluent
        annotation (Placement(transformation(extent={{88,-28},{108,-8}})));
      WasteWater.BSM1.SludgeSink WasteSludge
        annotation (Placement(transformation(extent={{87,-51},{107,-31}})));
      WasteWater.BSM1.divider2 divider
        annotation (Placement(transformation(extent={{20,-6},{40,14}})));
      WasteWater.BSM1.nitri tank3(V=1333)
        annotation (Placement(transformation(extent={{-6,-6},{14,14}})));
      WasteWater.BSM1.nitri tank2(V=1333)
        annotation (Placement(transformation(extent={{-34,-6},{-14,14}})));
      WasteWater.BSM1.deni tank1(V=3000)
        annotation (Placement(transformation(extent={{-65,-6},{-45,14}})));
      WasteWater.BSM1.mixer3 mixer
        annotation (Placement(transformation(extent={{-104,22},{-84,42}})));
      Modelica.Blocks.Sources.CombiTimeTable CombiTableTime(
        fileName=Modelica.Utilities.Files.loadResource("modelica://WasteWater/Resources/ASM1/Inf_dry.txt"),
        table=[0,0; 1,1],
        columns=integer(({16,3,4,5,6,7,8,9,10,11,12,13,14,15})),
        tableName="Inf_dry",
        tableOnFile=("Inf_dry") <> "NoName") annotation (Placement(transformation(
              extent={{-114,78},{-94,98}})));
      WasteWater.BSM1.WWSource WWSource
        annotation (Placement(transformation(extent={{-88,78},{-68,98}})));
      WasteWater.BSM1.blower blower1(Q_max=34574.2654508612)
        annotation (Placement(transformation(extent={{-33,-62},{-13,-42}})));
      WasteWater.BSM1.blower blower2(Q_max=34574.2654508612)
        annotation (Placement(transformation(extent={{-6,-62},{14,-42}})));
      WasteWater.BSM1.sensor_O2 sensor_O2
        annotation (Placement(transformation(extent={{0,24},{20,44}})));
      Modelica.Blocks.Math.Feedback Feedback annotation (Placement(transformation(
              extent={{62,40},{82,60}})));
      Modelica.Blocks.Continuous.PI PI1(T=0.001, k=500, initType=Modelica.Blocks.Types.Init.InitialState)
        annotation (Placement(transformation(extent={{88,40},{108,60}})));
      Modelica.Blocks.Sources.Constant Constant1 annotation (Placement(
            transformation(extent={{-67,-87},{-47,-67}})));
      WasteWater.BSM1.pump RecyclePump(Q_max=46115) annotation (Placement(
            transformation(
            origin={-84,-12},
            extent={{-10,-10},{10,10}},
            rotation=180)));
      WasteWater.BSM1.pump ReturnPump(Q_max=9223) annotation (Placement(
            transformation(
            origin={26,-26},
            extent={{-10,-10},{10,10}},
            rotation=180)));
      WasteWater.BSM1.pump WastePump(Q_max=193)
        annotation (Placement(transformation(extent={{59,-55},{79,-35}})));
      Modelica.Blocks.Sources.Constant Constant2 annotation (Placement(
            transformation(extent={{22,-68},{42,-48}})));
      Modelica.Blocks.Sources.Constant Temperature(k=15)
        annotation (Placement(transformation(extent={{-94,50},{-82,62}})));
      sensor_NH sensor_NH1 annotation (Placement(transformation(extent={{64,15},{
                80,31}})));
      WasteWater.ASM1.sensor_NO sensor_NO1 annotation (Placement(transformation(
              extent={{81,15},{97,31}})));
      WasteWater.ASM1.sensor_TKN sensor_TKN1 annotation (Placement(transformation(
              extent={{97,14},{113,30}})));
      WasteWater.ASM1.sensor_COD sensor_COD1 annotation (Placement(transformation(
              extent={{97,-5},{113,11}})));
      Modelica.Blocks.Sources.Step OxygenSetpoint(height=1.5)
        annotation (Placement(transformation(extent={{37,40},{57,60}})));
      WasteWater.ASM1.SecClar.Krebs.SecClarModKrebs Settler annotation (Placement(
            transformation(extent={{48,-4},{68,16}})));
      WasteWater.ASM1.sensor_TSS sensor_TSS1 annotation (Placement(transformation(
              extent={{32,14},{49,30}})));
    equation
      connect(tank3.Out, divider.In) annotation (Line(points={{14,4},{17,4},{17,
              4.3},{20,4.3}}));
      connect(mixer.Out, tank1.In) annotation (Line(points={{-84,31.6},{-77,
              31.6},{-77,4},{-65,4}}));
      connect(mixer.In1, WWSource.Out) annotation (Line(points={{-104,35.5},{
              -104,74},{-68,74},{-68,81},{-68.2,81}}));
      connect(CombiTableTime.y, WWSource.data)
        annotation (Line(points={{-93,88},{-87,88}}));
      connect(blower2.AirOut, tank3.AirIn) annotation (Line(points={{3,-42},{3,
              -24},{4,-24},{4,-5.8}}));
      connect(Feedback.y, PI1.u) annotation (Line(points={{81,50},{86,50}}));
      connect(PI1.y, blower2.u) annotation (Line(points={{109,50},{114,50},{114,
              -84},{18,-84},{18,-55},{13.8,-55},{13.8,-55}}));
      connect(divider.Out2, RecyclePump.In) annotation (Line(points={{40,2.5},{
              44,2.5},{44,-8.7},{-74,-8.7}}));
      connect(RecyclePump.Out, mixer.In3) annotation (Line(points={{-94,-14.8},
              {-104,-14.8},{-104,27.5}}));
      connect(ReturnPump.Out, mixer.In2) annotation (Line(points={{16,-28.8},{
              15.5,-28.8},{15.5,-30},{-112,-30},{-112,31.5},{-104,31.5}}));
      connect(sensor_O2.So, Feedback.u2)
        annotation (Line(points={{19.8,34},{72,34},{72,42}}));
      connect(Temperature.y, tank1.T)
        annotation (Line(points={{-81.4,56},{-71,56},{-71,8},{-65,8}}, color={0,0,
              255}));
      connect(Temperature.y, tank2.T)
        annotation (Line(points={{-81.4,56},{-39,56},{-39,8},{-34,8}}, color={0,0,
              255}));
      connect(Temperature.y, tank3.T) annotation (Line(points={{-81.4,56},{-39,
              56},{-39,14},{-5.9,14},{-5.9,8},{-6,8}},
                                                     color={0,0,255}));
      connect(OxygenSetpoint.y, Feedback.u1)
        annotation (Line(points={{58,50},{64,50}}, color={0,0,255}));
      connect(Constant1.y, blower1.u) annotation (Line(points={{-46,-77},{-7.2,
              -77},{-7.2,-55},{-13.2,-55}}, color={0,0,255}));
      connect(WastePump.Out, WasteSludge.In) annotation (Line(points={{79,-42.2},
              {81,-42.2},{81,-42},{83,-42},{83,-42},{87,-42}}));
      connect(WastePump.u, Constant2.y)
        annotation (Line(points={{60.1,-42.5},{46,-42.5},{46,-58},{43,-58}},
                                                                       color={0,0,
              255}));
      connect(tank2.Out, tank3.In) annotation (Line(points={{-14,4},{-6,4}}));
      connect(tank1.Out, tank2.In) annotation (Line(points={{-45,4},{-34,4}}));
      connect(blower1.AirOut, tank2.AirIn) annotation (Line(points={{-24,-42},{
              -24,-24},{-24,-5.8},{-24,-5.8}}));
      connect(Constant1.y, RecyclePump.u) annotation (Line(points={{-46,-77},{
              -39,-77},{-39,-14.5},{-75.1,-14.5}},
                                           color={0,0,255}));
      connect(Settler.Effluent, Effluent.In) annotation (Line(points={{68.2,
              11.7},{78,11.7},{78,-16},{88,-16}}));
      connect(Settler.Return, ReturnPump.In) annotation (Line(points={{55,-3.6},
              {55,-22.7},{36,-22.7}}));
      connect(WastePump.In, Settler.Waste) annotation (Line(points={{59,-48.3},
              {52,-48.3},{52,-31},{61,-31},{61,-3.6}}));
      connect(sensor_NH1.In, Settler.Effluent) annotation (Line(points={{72,15},
              {72,11.7},{68.2,11.7}}));
      connect(sensor_NO1.In, Settler.Effluent) annotation (Line(points={{89,15},
              {89,11.7},{68.2,11.7}}));
      connect(sensor_TKN1.In, Settler.Effluent) annotation (Line(points={{105,14},
              {105,11.7},{68.2,11.7}}));
      connect(sensor_COD1.In, Settler.Effluent) annotation (Line(points={{105,-5},
              {105,11.7},{68.2,11.7}}));
      connect(Constant2.y, ReturnPump.u) annotation (Line(points={{43,-58},{48,
              -58},{48,-28.5},{34.9,-28.5}}, color={0,0,255}));
      connect(divider.Out1, Settler.Feed) annotation (Line(points={{40,6.6},{
              44.5,6.6},{44.5,7.4},{48,7.4}}));
      connect(tank3.MeasurePort, sensor_O2.In) annotation (Line(points={{9.5,8.5},
              {9.5,16},{10,16},{10,25}}));
      connect(sensor_TSS1.In, divider.Out1) annotation (Line(points={{40.5,14},
              {40.5,10.6},{40,10.6},{40,6}}));

      annotation (
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-120,-100},{120,105}},
            grid={1,1}), graphics={Line(points={{-22,58},{-22,58}}, color={0,0,
                  255})}),
        Documentation(info="This fictitious plant provides an ASM1 example model with a small number of equations.
It consists of one denitrification and 2 nitrification tanks and a settler.

Change into the directory ../ASM1 and translate the model.
Before simulating the model load initial values from the script file small_asm1.mos
that is provided besides the model.
A 14 days dynamic influent data file is provided. So you may simulate up to 14 days.
But start with 1 day as it may take some time for simulation.
After simulation you may have a look at internal concentrations but most interesting
are the relevant concentrations at the effluent of a plant which can be viewed via the
sensors at the effluent of the secondary clarifier.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de
"));
    end SmallPlant;

    class BenchPlant "COST Benchmark WWTP Configuration"
      import WasteWater;

      //Q_air=34574.2654508612 is equal to a Kla of 10 h^-1 from COST benchmark
      //Q_air=12100.99290780142 is equal to a Kla of 3.5 h^-1 from COST benchmark
      extends Modelica.Icons.Example;

      WasteWater.BSM1.EffluentSink Effluent
        annotation (Placement(transformation(extent={{79,-25.5},{99,-5.5}})));
      WasteWater.BSM1.SludgeSink WasteSludge
        annotation (Placement(transformation(extent={{120,-69},{140,-49}})));
      WasteWater.BSM1.SecClarModTakacs Settler
        annotation (Placement(transformation(extent={{33,-5},{53,14}})));
      WasteWater.BSM1.divider2 divider
        annotation (Placement(transformation(extent={{8,-6},{28,14}})));
      WasteWater.BSM1.nitri5
                 tank5(V=1333) annotation (Placement(transformation(extent={{-18,-6},
                {2,14}})));
      WasteWater.BSM1.nitri tank4(V=1333)
        annotation (Placement(transformation(extent={{-44,-6},{-24,14}})));
      WasteWater.BSM1.nitri tank3(V=1333)
        annotation (Placement(transformation(extent={{-70,-6},{-50,14}})));
      WasteWater.BSM1.deni tank2
        annotation (Placement(transformation(extent={{-65.5,22},{-45.5,42}})));
      WasteWater.BSM1.deni tank1
        annotation (Placement(transformation(extent={{-93,22},{-73,42}})));
      WasteWater.BSM1.mixer3 mixer
        annotation (Placement(transformation(extent={{-120,22.5},{-100,42.5}})));
      Modelica.Blocks.Sources.CombiTimeTable CombiTableTime(
        fileName=Modelica.Utilities.Files.loadResource("modelica://WasteWater/Resources/ASM1/Inf_dry.txt"),
        table=[0,0; 1,1],
        columns=integer(({2,3,4,5,6,7,8,9,10,11,12,13,14,15})),
        tableName="Inf_dry",
        tableOnFile=("Inf_dry") <> "NoName") annotation (Placement(transformation(
              extent={{-80,60},{-100,80}})));
      WasteWater.BSM1.WWSource WWSource
        annotation (Placement(transformation(extent={{-110,60},{-130,80}})));
      WasteWater.BSM1.sensor_NO sensor_NO
        annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
      WasteWater.BSM1.sensor_O2 sensor_O2
        annotation (Placement(transformation(extent={{-11,13},{6,30}})));
      Modelica.Blocks.Sources.Constant Constant1(k=55338)
                                                 annotation (Placement(
            transformation(extent={{-5,-5},{5,5}},
            rotation=180,
            origin={5,-20})));
      WasteWater.BSM1.pump RecyclePump(Q_max=55338) annotation (Placement(
            transformation(
            origin={-30,-32.5},
            extent={{-10,10.5},{10,-10.5}},
            rotation=180)));
      WasteWater.BSM1.pump ReturnPump(Q_max=18446) annotation (Placement(
            transformation(
            origin={-30,-47},
            extent={{-10,-10},{10,10}},
            rotation=180)));
      WasteWater.BSM1.pump WastePump(Q_max=385)
        annotation (Placement(transformation(extent={{56,-37},{76,-57}})));
      Modelica.Blocks.Sources.Constant Constant2 annotation (Placement(
            transformation(extent={{-5,-5},{5,5}},
            rotation=180,
            origin={5,-60})));
      Modelica.Blocks.Sources.Constant Temperature(k=15)
        annotation (Placement(transformation(extent={{-106,45},{-96,55}})));
      sensor_NH sensor_NH1 annotation (Placement(transformation(extent={{45,14},
                {61,30}})));
      WasteWater.ASM1.sensor_NO sensor_NO1 annotation (Placement(transformation(
              extent={{60,14},{76,30}})));
      WasteWater.ASM1.sensor_TKN sensor_TKN1 annotation (Placement(transformation(
              extent={{75,14},{91,30}})));
      WasteWater.ASM1.sensor_COD sensor_COD1 annotation (Placement(transformation(
              extent={{90,14},{106,30}})));
      WasteWater.ASM1.sensor_TSS sensor_TSS1 annotation (Placement(transformation(
              extent={{20,15},{36,30}})));
      WasteWater.BSM1.sensor_effluent_quality sensor_effluent_quality
        annotation (Placement(transformation(extent={{57,-10},{77,10}})));
      WasteWater.BSM1.sensor_pump_energy sensor_pump_energy
        annotation (Placement(transformation(extent={{30,-40},{50,-20}})));
      WasteWater.BSM1.sensor_influent_quality sensor_influent_quality1
        annotation (Placement(transformation(extent={{-130,40},{-110,60}})));
      WasteWater.BSM1.sensor_sludge_production sensor_sludge_production
        annotation (Placement(transformation(extent={{85,-60},{105,-40}})));
      WasteWater.BSM1.OCI oCI
        annotation (Placement(transformation(extent={{129,-38},{149,-18}})));
      WasteWater.BSM1.sensor_aeration_energy sensor_aeration_energy
        annotation (Placement(transformation(extent={{-20.5,40},{-0.5,60}})));
    equation
      connect(divider.Out1, Settler.Feed) annotation (Line(points={{28,6.5},{33,
              6.5},{33,5.83}}));
      connect(tank5.Out, divider.In) annotation (Line(points={{2,4},{7,4},{7,
              4.1},{8,4.1}}));
      connect(tank4.Out, tank5.In) annotation (Line(points={{-24,4},{-18,4}}));
      connect(tank3.Out, tank4.In) annotation (Line(points={{-50,4},{-44,4}}));
      connect(tank3.In, tank2.Out) annotation (Line(points={{-70,4},{-80,4},{
              -80,20},{-40,20},{-40,32},{-45.5,32}}));
      connect(tank1.Out, tank2.In) annotation (Line(points={{-73,32},{-65.5,32}}));
      connect(mixer.Out, tank1.In) annotation (Line(points={{-100,32.1},{-97,
              32.1},{-97,32},{-93,32}}));
      connect(sensor_NO.In, tank2.MeasurePort) annotation (Line(points={{-50,40},
              {-50,38},{-50,36.5}}));
      connect(RecyclePump.Out, mixer.In3) annotation (Line(points={{-40,-29.56},
              {-120,-29.56},{-120,28}}));
      connect(ReturnPump.Out, mixer.In2) annotation (Line(points={{-40,-49.8},{
              -40,-49.8},{-130,-49.8},{-130,32},{-120,32}}));
      connect(Temperature.y, tank1.T)
        annotation (Line(points={{-95.5,50},{-93.5,50},{-93.5,36},{-93,36}},
                                                                           color=
              {0,0,255}));
      connect(Temperature.y, tank2.T)
        annotation (Line(points={{-95.5,50},{-70,50},{-70,36},{-65.5,36}},
                                                                         color={0,
              0,255}));
      connect(Temperature.y, tank4.T) annotation (Line(points={{-95.5,50},{-70,
              50},{-70,15},{-44,15},{-44,8}},
                                          color={0,0,255}));
      connect(sensor_NH1.In, Settler.Effluent) annotation (Line(points={{53,14},
              {53,9.915},{53.2,9.915}}));
      connect(sensor_NO1.In, Settler.Effluent) annotation (Line(points={{68,14},
              {68,9.915},{53.2,9.915}}));
      connect(sensor_TKN1.In, Settler.Effluent) annotation (Line(points={{83,14},
              {83,9.915},{53.2,9.915}}));
      connect(sensor_COD1.In, Settler.Effluent) annotation (Line(points={{98,14},
              {98,9.915},{53.2,9.915}}));
      connect(tank5.MeasurePort, sensor_O2.In) annotation (Line(points={{-2.5,
              8.5},{-2.5,13}}));
      connect(sensor_TSS1.In, divider.Out1) annotation (Line(points={{28,15},{
              28,11},{28,6.5}}));

      connect(CombiTableTime.y, WWSource.data) annotation (Line(points={{-101,70},
              {-111,70}},                  color={0,0,127}));
      connect(sensor_effluent_quality.In, Settler.Effluent) annotation (Line(
          points={{60,0},{60,9.915},{53.2,9.915}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(Effluent.In, sensor_effluent_quality.Out) annotation (Line(
          points={{79,-13.5},{67,-13.5},{67,-10}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(divider.Out2, sensor_pump_energy.In_a) annotation (Line(
          points={{28,2.4},{28,-5},{34,-5},{34,-24}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(Settler.Return, sensor_pump_energy.In_r) annotation (Line(
          points={{40,-4.62},{40,-24}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(Settler.Waste, sensor_pump_energy.In_w) annotation (Line(
          points={{46,-4.62},{46,-24}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(sensor_pump_energy.Out_w, WastePump.In) annotation (Line(
          points={{46,-36},{46,-43.7},{56,-43.7}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(WWSource.Out, sensor_influent_quality1.In) annotation (Line(
          points={{-129.8,63},{-130,63},{-130,50},{-127,50}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(sensor_influent_quality1.Out, mixer.In1) annotation (Line(
          points={{-120,41},{-120,36}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(sensor_sludge_production.Out, WasteSludge.In) annotation (Line(
          points={{95,-60},{95,-60.2},{120,-60.2}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(sensor_sludge_production.SP, oCI.SP) annotation (Line(
          points={{104.8,-50},{110,-50},{110,-34},{131,-34}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(sensor_pump_energy.PE, oCI.PE) annotation (Line(
          points={{49.8,-30},{131,-30}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(sensor_effluent_quality.EQ, oCI.EC) annotation (Line(
          points={{76.8,0},{110,0},{110,-26},{131,-26}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(tank5.Kla, sensor_aeration_energy.Kla5) annotation (Line(
          points={{-11,9},{-11,22.65},{-11.1,22.65},{-11.1,42}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(tank5.T, tank4.T) annotation (Line(
          points={{-18,8},{-18,15},{-44,15},{-44,8}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(tank3.T, tank4.T) annotation (Line(
          points={{-70,8},{-70,15},{-44,15},{-44,8}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(sensor_sludge_production.In, WastePump.Out) annotation (Line(
          points={{88,-50},{77,-50},{77,-49.8},{76,-49.8}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(Constant2.y, WastePump.u) annotation (Line(
          points={{-0.5,-60},{-10,-60},{-10,-49.5},{57.1,-49.5}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(ReturnPump.u, WastePump.u) annotation (Line(
          points={{-21.1,-49.5},{57.1,-49.5}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(sensor_pump_energy.Out_a, RecyclePump.In) annotation (Line(
          points={{34,-36},{9,-36},{9,-35.965},{-20,-35.965}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(RecyclePump.u, Constant1.y) annotation (Line(
          points={{-21.1,-29.875},{-10.05,-29.875},{-10.05,-20},{-0.5,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(oCI.AE, sensor_aeration_energy.AE) annotation (Line(
          points={{131,-22},{120,-22},{120,50},{-1.5,50}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(sensor_pump_energy.Out_r, ReturnPump.In) annotation (Line(
          points={{40,-36},{40,-43.7},{-20,-43.7}},
          color={0,0,0},
          smooth=Smooth.None));
      annotation (
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-140,-100},{150,105}},
            grid={1,1}), graphics),
        Documentation(info="This ASM1 plant consists of 2 denitrification tanks (tank1 and tank2),
3 nitrification tanks (tank3 - tank5) and a secondary clarifier by Takacs.
Furthermore there are 2 control loops modelled.
This configuration corresponds to the COST simulation benchmark [1].

Change into the directory ../ASM1 and translate the model.
Before simulating the model load initial values from the script file bench_asm1.mos
that is provided besides the model.
A 14 days dynamic influent data file is provided. So you may simulate up to 14 days.
But start with 1 day as it may take some time for simulation.
After simulation you may have a look at internal concentrations but most interesting
are the relevant concentrations at the effluent of a plant which can be viewed via the
sensors at the effluent of the secondary clarifier.

References:

[1] J.B. Copp: The COST Simulation Benchmark. 2000. http://www.ensic.u-nancy.fr/COSTWWTP/


PS: For those who want to reproduce the exact figures from the COST simulation benchmark some remarks:
    The aeration system in this library is different from that in COST, so be sure to produce an airflow
    corresponding to the desired Kla in COST. Furthermore in this library biological parameters are standard
    parameters from the ASM1 distribution and implemented with temperature dependency which may vary a bit from
    the parameter set used in COST.
    But it is possible. During the validation phase of this library the steady state and dynamic results
    from the COST simulation benchmark could exactly be reproduced.
"),
        experiment(
          StopTime=14,
          __Dymola_NumberOfIntervals=10000,
          __Dymola_Algorithm="lsodar"),
        __Dymola_experimentSetupOutput);
    end BenchPlant;

    model ComplexPlant "Complex ASM1 WWTP"
      import WasteWater;

      extends Modelica.Icons.Example;
      ControlledDivider2 cdivider1 annotation (Placement(transformation(extent={{
                -168,65},{-148,85}})));
      Modelica.Blocks.Sources.Constant Constant2(k=0.8)
        annotation (Placement(transformation(extent={{-178,52},{-168,62}})));
      blower blower1(Q_max=162816)
        "there exist 4 blowers of 4240 Nm3/h each, Q_max adusted according active aerated tanks"
         annotation (Placement(transformation(extent={{145,-16},{165,4}})));
      nitri nitri2(
        V=2772,
        alpha=0.305,
        de=5.24,
        R_air=20) annotation (Placement(transformation(extent={{110,18},{130,38}})));
      deni anaerob(V=1287) annotation (Placement(transformation(extent={{-138,13},
                {-118,33}})));
      deni deni1(V=2772) annotation (Placement(transformation(extent={{-80,14},{
                -60,34}})));
      deni deni3(V=2772) annotation (Placement(transformation(extent={{80,18},{
                100,38}})));
      deni deni2(V=2772) annotation (Placement(transformation(extent={{-20,14},{0,
                34}})));
      nitri nitri3(
        V=5602,
        alpha=0.305,
        de=5.24,
        R_air=21) annotation (Placement(transformation(extent={{144,18},{164,38}})));
      blower blower2(Q_max=81408)
        "there exist 4 blowers of 4240 Nm3/h each, Q_max adjusted according active aerated tanks"
         annotation (Placement(transformation(extent={{111,-13},{131,7}})));
      pump ReturnPump(Q_max=60480) annotation (Placement(transformation(
            origin={-44,-94},
            extent={{-10,-10},{10,10}},
            rotation=180)));
      pump RecyclePump(Q_max=60480) annotation (Placement(transformation(
            origin={10,-37},
            extent={{-10,-10},{10,10}},
            rotation=180)));
      pump WastePump(Q_max=1920) annotation (Placement(transformation(extent={{
                128,-104},{148,-84}})));
      ControlledDivider2 cdivider2 annotation (Placement(transformation(extent={{
                -42,68},{-22,88}})));
      EffluentSink Effluent annotation (Placement(transformation(extent={{170,-72},
                {190,-52}})));
      SludgeSink WasteSludge annotation (Placement(transformation(extent={{170,
                -100},{190,-80}})));
      mixer2 mixer2_1 annotation (Placement(transformation(extent={{-45,16},{-25,
                36}})));
      mixer2 mixer2_2 annotation (Placement(transformation(extent={{50,17},{70,37}})));
      mixer3 mixer3_1 annotation (Placement(transformation(extent={{-107,14},{-87,
                34}})));
      mixer2 mixer2_5 annotation (Placement(transformation(extent={{-140,-15},{
                -120,5}})));
      divider2 divider2_1 annotation (Placement(transformation(
            origin={66,-32},
            extent={{-10,-10},{10,10}},
            rotation=180)));
      ControlledDivider2 cdivider3 annotation (Placement(transformation(
            origin={-122,-87},
            extent={{-10,-10},{10,10}},
            rotation=180)));
      Modelica.Blocks.Sources.Constant Constant4(k=0.5)
        annotation (Placement(transformation(extent={{-148,-80},{-138,-70}})));
      nitri nitri1(
        V=5602,
        alpha=0.305,
        de=5.24,
        R_air=21) annotation (Placement(transformation(extent={{12,14},{32,34}})));
      blower blower3(Q_max=162816)
        "there exist 4 blowers of max 4240 Nm3/h, Q_max adusted according active aerated tanks"
         annotation (Placement(transformation(extent={{13,-23},{33,-3}})));
      Modelica.Blocks.Sources.Constant Constant7(k=0.56)
        annotation (Placement(transformation(extent={{-66,98},{-56,108}})));
      mixer2 mixer2_3 annotation (Placement(transformation(extent={{-168,18},{
                -148,38}})));
      PreClar.preclar3 Preclaryfier(V=1372, n_corr=2.138) annotation (Placement(
            transformation(extent={{-136,68},{-116,88}})));
      FlowSource FlowInput annotation (Placement(transformation(extent={{-176,94},
                {-156,114}})));
      ControlledDivider2 ControlledDivider2_1 annotation (Placement(
            transformation(
            origin={-41,-39},
            extent={{-10,-10},{10,10}},
            rotation=180)));
      Modelica.Blocks.Sources.Constant Constant6 annotation (Placement(
            transformation(extent={{-46,-10},{-36,0}})));
      sensor_NO sensor_NO1 annotation (Placement(transformation(extent={{-14,40},
                {6,60}})));
      sensor_NO sensor_NO3 annotation (Placement(transformation(extent={{90,-56},
                {110,-36}})));
      sensor_NH sensor_NH2 annotation (Placement(transformation(extent={{112,-56},
                {132,-36}})));
      sensor_TSS sensor_TSS1 annotation (Placement(transformation(extent={{-6,-86},
                {14,-66}})));
      Modelica.Blocks.Sources.Constant Temperature(k=11.5)
        annotation (Placement(transformation(extent={{58,44},{78,64}})));
      WasteWater.Misc.RecycleController2 RecycleController1(NO3min=1.5)
        annotation (Placement(transformation(extent={{15,-63},{29,-49}})));
      WasteWater.Misc.ReturnController ReturnController1 annotation (Placement(
            transformation(extent={{-28,-115},{-14,-101}})));
      sensor_Q sensor_Q1 annotation (Placement(transformation(extent={{-194,66},{
                -174,86}})));
      WasteWater.Misc.TwoPoint TwoPoint1(
        on=4.5,
        off=4.0,
        out_on=2.5,
        out_off=1.5) annotation (Placement(transformation(extent={{-8,100},{2,110}})));
      WasteWater.Misc.TwoPoint TwoPoint2(
        on=4.5,
        off=4.0,
        out_on=2.5,
        out_off=1.5) annotation (Placement(transformation(extent={{110,100},{120,
                110}})));
      WasteWater.ASM1.sensor_NH sensor_NH1 annotation (Placement(transformation(
              extent={{142,38},{162,58}})));
      Modelica.Blocks.Math.Feedback Feedback1 annotation (Placement(
            transformation(extent={{12,100},{22,110}})));
      Modelica.Blocks.Math.Feedback Feedback2 annotation (Placement(
            transformation(extent={{127,99},{137,109}})));
      WasteWater.ASM1.sensor_O2 sensor_O2_1 annotation (Placement(transformation(
              extent={{164,38},{184,58}})));
      WasteWater.ASM1.sensor_O2 sensor_O2_2 annotation (Placement(transformation(
              extent={{14,38},{34,58}})));
      WasteWater.ASM1.sensor_COD sensor_COD1 annotation (Placement(transformation(
              extent={{134,-56},{154,-36}})));
      WasteWater.ASM1.sensor_COD sensor_COD2 annotation (Placement(transformation(
              extent={{-118,94},{-98,114}})));
      WasteWater.ASM1.Examples.JenaSecClarModTakacs Settler(hsc=3.46, Asc=3704)
        "The depth is calculated based on V and A of the settler and not the true depth."
         annotation (Placement(transformation(extent={{68,-84},{88,-64}})));
      WasteWater.ASM1.sensor_TKN sensor_TKN1 annotation (Placement(transformation(
              extent={{68,-56},{88,-36}})));
      Modelica.Blocks.Sources.CombiTimeTable CombiTableTime1(
        fileName=Modelica.Utilities.Files.loadResource("modelica://WasteWater/Resources/ASM1/drysim130303.txt"),
        table=[0,0; 1,1],
        columns=integer(({2})),
        tableName="drysim130303",
        tableOnFile=("drysim130303") <> "NoName")
        annotation (Placement(transformation(extent={{-198,96},{-184,112}})));
      Modelica.Blocks.Sources.CombiTimeTable CombiTableTime2(
        fileName=Modelica.Utilities.Files.loadResource("modelica://WasteWater/Resources/ASM1/drysim130303.txt"),
        table=[0,0; 1,1],
        columns=integer(({3,7,4,6})),
        tableName="drysim130303",
        tableOnFile=("drysim130303") <> "NoName")
        annotation (Placement(transformation(extent={{-134,96},{-118,112}})));
      Modelica.Blocks.Nonlinear.FixedDelay FixedDelay1(delayTime=1/24/6)
        annotation (Placement(transformation(extent={{-16,-63},{-2,-49}})));
      Modelica.Blocks.Math.Feedback Feedback3 annotation (Placement(
            transformation(
            origin={84,-13},
            extent={{5,-5},{-5,5}},
            rotation=180)));
      WasteWater.ASM1.sensor_O2 sensor_O2_3 annotation (Placement(transformation(
              extent={{116,40},{134,58}})));
      Modelica.Blocks.Math.Gain Gain1(k=500) annotation (Placement(transformation(
              extent={{29,95},{49,115}})));
      Modelica.Blocks.Math.Gain Gain2(k=500) annotation (Placement(transformation(
              extent={{144,93},{164,113}})));
      Modelica.Blocks.Math.Gain Gain3(k=500) annotation (Placement(transformation(
              extent={{94,-28},{108,-14}})));
      WasteWater.Misc.TwoPoint TwoPoint3(
        on=4.5,
        off=4.0,
        out_on=2.0,
        out_off=1.0) annotation (Placement(transformation(extent={{63,-19},{75,-7}})));
      Modelica.Blocks.Sources.Step Step1(
        height=0.125,
        offset=-0.5,
        startTime=2.375) annotation (Placement(transformation(extent={{100,-90},{
                110,-80}})));
      WasteWater.ASM1.sensor_TSS sensor_TSS2 annotation (Placement(transformation(
              extent={{166,5},{184,24}})));
    equation
      connect(deni3.Out, nitri2.In) annotation (Line(points={{100,28},{110,28}}));
      connect(nitri2.Out, nitri3.In) annotation (Line(points={{130,28},{144,28}}));
      connect(nitri3.Out, divider2_1.In) annotation (Line(points={{164,28},{190,
              28},{190,-32.3},{76,-32.3}}));
      connect(anaerob.Out, mixer3_1.In2) annotation (Line(points={{-118,23},{
              -110.5,23},{-110.5,24.5},{-107,24.5}}));
      connect(mixer2_5.Out, mixer3_1.In3) annotation (Line(points={{-120,-4.6},{
              -111,-4.6},{-111,20.5},{-107,20.5}}));
      connect(Constant4.y, cdivider3.u)
        annotation (Line(points={{-137.5,-75},{-122,-75},{-122,-81}}, color={0,0,
              255}));
      connect(nitri1.Out, mixer2_2.In2) annotation (Line(points={{33,24},{39,24},
              {39,25.5},{51,25.5}}));
      connect(deni2.Out, nitri1.In) annotation (Line(points={{5.55112e-16,24},{12,
              24}}));
      connect(cdivider2.Out1, mixer2_2.In1) annotation (Line(points={{-22,80.5},{
              40,80.5},{40,29.5},{50,29.5}}));
      connect(cdivider3.Out2, mixer2_5.In1) annotation (Line(points={{-132,-85.5},
              {-172,-85.5},{-172,-2.5},{-140,-2.5}}));
      connect(deni1.Out, mixer2_1.In1) annotation (Line(points={{-60,24},{-56,24},
              {-56,29},{-50.75,29},{-50.75,28.5},{-45,28.5}}));
      connect(Constant7.y, cdivider2.u) annotation (Line(points={{-55.5,103},{
              -50.5,103},{-50.5,54},{-32,54},{-32,72}}, color={0,0,255}));
      connect(cdivider3.In, ReturnPump.Out) annotation (Line(points={{-112,-87.3},
              {-86.1,-87.3},{-86.1,-96.8},{-54.2,-96.8}}));
      connect(cdivider3.Out1, mixer2_3.In2) annotation (Line(points={{-132,-89.5},
              {-184,-89.5},{-184,26.5},{-168,26.5}}));
      connect(cdivider1.Out2, mixer2_3.In1) annotation (Line(points={{-147,73.5},
              {-143,73.5},{-143,40},{-184,40},{-184,30.5},{-168,30.5}}));
      connect(WastePump.Out, WasteSludge.In) annotation (Line(points={{148.2,
              -91.2},{170,-91.2}}));
      connect(ControlledDivider2_1.Out2, mixer2_1.In2) annotation (Line(points={{
              -51,-38.5},{-52,-38.5},{-52,24.5},{-45,24.5}}));
      connect(ControlledDivider2_1.Out1, mixer2_5.In2) annotation (Line(points={{
              -51,-42.5},{-150,-42.5},{-150,-6.5},{-140,-6.5}}));
      connect(Constant6.y, ControlledDivider2_1.u) annotation (Line(points={{-36,
              -5},{-31,-5},{-31,-18},{-41,-18},{-41,-34}}, color={0,0,255}));
      connect(Preclaryfier.In, cdivider1.Out1) annotation (Line(points={{-136,78},
              {-146.5,78},{-146.5,77.5},{-147,77.5}}));
      connect(mixer2_2.Out, deni3.In) annotation (Line(points={{70,27.4},{70,28},
              {80,28}}));
      connect(sensor_NO1.In, deni2.MeasurePort) annotation (Line(points={{-4,40},
              {-4,28.5},{-4.5,28.5}}));
      connect(Temperature.y, anaerob.T) annotation (Line(points={{80,54},{94,54},
              {94,34},{-138,34},{-138,27}}, color={0,0,255}));
      connect(Temperature.y, deni1.T) annotation (Line(points={{80,54},{94,54},{
              94,34},{-80,34},{-80,28}}, color={0,0,255}));
      connect(Temperature.y, deni2.T) annotation (Line(points={{80,54},{94,54},{
              94,34},{-20,34},{-20,28}}, color={0,0,255}));
      connect(Temperature.y, nitri1.T) annotation (Line(points={{80,54},{94,54},{
              94,34},{12,34},{12,28},{12.2,28}}, color={0,0,255}));
      connect(Temperature.y, deni3.T)
        annotation (Line(points={{80,54},{94,54},{94,34},{80,34},{80,32}}, color=
              {0,0,255}));
      connect(Temperature.y, nitri2.T) annotation (Line(points={{80,54},{94,54},{
              94,34},{110.5,34},{110.5,32},{110,32}}, color={0,0,255}));
      connect(Temperature.y, nitri3.T) annotation (Line(points={{80,54},{94,54},{
              94,34},{144,34},{144,32}}, color={0,0,255}));
      connect(RecycleController1.out, RecyclePump.u) annotation (Line(points={{
              29.7,-56},{34,-56},{34,-38.5},{18.8,-38.5}}, color={0,0,255}));
      connect(FlowInput.Out, sensor_Q1.In) annotation (Line(points={{-156,97},{
              -148,97},{-148,90},{-194,90},{-194,76}}));
      connect(sensor_Q1.Q, ReturnController1.in1) annotation (Line(points={{-184,
              66},{-184,65.5},{-194,65.5},{-194,-108},{-28.7,-108}}, color={0,0,
              255}));
      connect(sensor_NH1.In, nitri3.MeasurePort) annotation (Line(points={{152,38},
              {152,34},{160,34},{160,32.5}}));
      connect(sensor_O2_2.In, nitri1.MeasurePort) annotation (Line(points={{24,38},
              {24,28.5},{28,28.5}}));
      connect(Feedback1.u2, sensor_O2_2.So)
        annotation (Line(points={{17,101},{17,83.5},{34,83.5},{34,48}}, color={0,
              0,255}));
      connect(sensor_O2_1.In, nitri3.MeasurePort) annotation (Line(points={{174,
              38},{174,34},{160.25,34},{160.25,32},{160,32}}));
      connect(Feedback2.u2, sensor_O2_1.So)
        annotation (Line(points={{132,101},{132,81},{184,81},{184,48}}, color={0,
              0,255}));
      connect(sensor_NH1.Snh, TwoPoint1.e) annotation (Line(points={{162,48},{162,
              72},{-14,72},{-14,105},{-8,105}}, color={0,0,255}));
      connect(TwoPoint2.e, sensor_NH1.Snh) annotation (Line(points={{110,105},{
              100,105},{100,72},{162,72},{162,48}}, color={0,0,255}));
      connect(sensor_COD2.In, Preclaryfier.In) annotation (Line(points={{-108,94},
              {-108,90},{-136,90},{-136,78}}));
      connect(Settler.Effluent, Effluent.In) annotation (Line(points={{88.2,-68.3},
              {108.65,-68.3},{108.65,-68},{128,-68},{128,-60},{170,-60}}));
      connect(sensor_NO3.In, Settler.Effluent) annotation (Line(points={{100,-56},
              {100,-68.3},{88.2,-68.3}}));
      connect(sensor_NH2.In, Settler.Effluent) annotation (Line(points={{122,-56},
              {122,-68.3},{88.2,-68.3}}));
      connect(sensor_COD1.In, Settler.Effluent) annotation (Line(points={{144,-56},
              {144,-60},{128,-60},{128,-68.3},{88.2,-68.3}}));
      connect(WastePump.In, Settler.Waste) annotation (Line(points={{128,-97.3},{
              82,-97.3},{82,-84}}));
      connect(Settler.Return, ReturnPump.In) annotation (Line(points={{74,-84},{
              74,-90.7},{-34,-90.7}}));
      connect(sensor_TSS1.In, Settler.Return) annotation (Line(points={{4,-86},{4,
              -91},{74,-91},{74,-84}}));
      connect(sensor_TKN1.In, Settler.Effluent) annotation (Line(points={{78,-56},
              {88,-56},{88,-68}}));
      connect(CombiTableTime1.y[1], FlowInput.data)
        annotation (Line(points={{-183.3,104},{-176,104}}, color={0,0,255}));
      connect(Preclaryfier.MeasurePort, CombiTableTime2.y) annotation (Line(
            points={{-122,88},{-122,94},{-116,94},{-116,104},{-117.2,104}}, color=
             {0,0,255}));
      connect(FixedDelay1.u, sensor_NO1.Sno) annotation (Line(points={{-18,-56},{
              -22,-56},{-22,-8},{10,-8},{10,50},{6,50}}, color={0,0,255}));
      connect(sensor_O2_3.In, nitri2.MeasurePort) annotation (Line(points={{125,
              40},{125,32.5},{126,32.5}}));
      connect(sensor_Q1.Out, cdivider1.In) annotation (Line(points={{-174,76},{
              -168,76}}));
      connect(cdivider2.Out2, mixer3_1.In1) annotation (Line(points={{-22,76},{
              -18,76},{-18,47},{-111,47},{-111,28.5},{-107,28.5}}));
      connect(mixer3_1.Out, deni1.In) annotation (Line(points={{-87,24},{-81,24}}));
      connect(mixer2_1.Out, deni2.In) annotation (Line(points={{-24,27},{-22.5,27},
              {-22.5,23},{-21,23}}));
      connect(FixedDelay1.y, RecycleController1.in1)
        annotation (Line(points={{-1,-56},{14.3,-56}}, color={0,0,255}));
      connect(Feedback2.u1, TwoPoint2.u)
        annotation (Line(points={{128,105},{120,105}}, color={0,0,255}));
      connect(TwoPoint1.u, Feedback1.u1)
        annotation (Line(points={{2.5,105},{13,105}}, color={0,0,255}));
      connect(mixer2_3.Out, anaerob.In) annotation (Line(points={{-148,29},{-143,
              29},{-143,23},{-138,23}}));
      connect(Gain1.u, Feedback1.y)
        annotation (Line(points={{27,104},{24,104},{24,105},{21,105}}, color={0,0,
              255}));
      connect(Gain2.u, Feedback2.y)
        annotation (Line(points={{141,103},{139,103},{139,104},{137,104}}, color=
              {0,0,255}));
      connect(Feedback3.y, Gain3.u) annotation (Line(points={{88.5,-13},{90.55,
              -13},{90.55,-21},{92.6,-21}}, color={0,0,255}));
      connect(sensor_O2_3.So, Feedback3.u2) annotation (Line(points={{134,49},{
              139,49},{139,16},{83,16},{83,-9}}, color={0,0,255}));
      connect(TwoPoint3.u, Feedback3.u1)
        annotation (Line(points={{75,-13},{80,-13}}, color={0,0,255}));
      connect(TwoPoint3.e, sensor_NH1.Snh) annotation (Line(points={{63,-13},{56,
              -13},{56,6},{104,6},{104,72},{162,72},{162,49}}, color={0,0,255}));
      connect(Step1.y, WastePump.u) annotation (Line(points={{111,-85},{119,-85},
              {119,-91},{129.2,-91}}, color={0,0,255}));
      connect(blower1.AirOut, nitri3.AirIn) annotation (Line(points={{155,4},{155,
              5.5},{154,5.5},{154,18}}));
      connect(Constant2.y, cdivider1.u)
        annotation (Line(points={{-167.5,57},{-158,57},{-158,69}}, color={0,0,255}));
      connect(Gain1.y, blower3.u) annotation (Line(points={{50,105},{55,105},{55,
              83},{49,83},{49,-15},{32,-15}}, color={0,0,255}));
      connect(Gain3.y, blower2.u) annotation (Line(points={{108.7,-21},{139,-21},
              {139,-6},{130.8,-6}}, color={0,0,255}));
      connect(blower1.u, Gain2.y)
        annotation (Line(points={{165.8,-9},{195,-9},{195,103},{165,103}}, color=
              {0,0,255}));
      connect(Settler.Feed, divider2_1.Out1) annotation (Line(points={{68,-72},{
              48,-72},{48,-35},{55,-35}}));
      connect(ReturnPump.u, ReturnController1.out) annotation (Line(points={{
              -35.1,-97},{-5,-97},{-5,-108},{-13,-108}}, color={0,0,255}));
      connect(RecyclePump.Out, ControlledDivider2_1.In) annotation (Line(points={
              {-8.88178e-16,-40},{-30,-40}}));
      connect(RecyclePump.In, divider2_1.Out2) annotation (Line(points={{20,-33.7},
              {38,-33.7},{38,-31},{56,-31}}));
      connect(blower3.AirOut, nitri1.AirIn) annotation (Line(points={{22,-3},{22,
              14}}));
      connect(cdivider2.In, Preclaryfier.Out) annotation (Line(points={{-43,78},{
              -115,78}}));
      connect(blower2.AirOut, nitri2.AirIn) annotation (Line(points={{120,7},{120,
              18}}));
      connect(sensor_TSS2.In, nitri3.Out) annotation (Line(points={{175,5},{175,
              28},{163,28}}));

      annotation (
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-200,-120},{200,120}},
            grid={1,1}), graphics),
        Documentation(info="This ASM1 example plant configuration is from a real municipal wastewater treatment plant
with a size of 145.000 p.e. It is a cascade-type continuous flow plant for a mean dry
weather inflow of 28.500 m3/d. It consists of a preclarifier, an anaerobic tank,
3 denitrification and 3 nitrification tanks and a secondary settler.
This model is an example for the Wastewater library and is not adapted with its parameters
to the reality, therefore simulation results do not reflect the real plant behaviour.

Change into the directory ../ASM1 and translate the model.
Before simulating the model load initial values from the script file complex_asm1.mos
that is provided besides the model.
A 14 days dynamic influent data file is provided. So you may simulate up to 14 days.
But start with 1 day as it may take some time for simulation.
After simulation you may have a look at internal concentrations but most interesting
are the relevant concentrations at the effluent of a plant which can be viewed via the
sensors at the effluent of the secondary clarifier.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de
"));
    end ComplexPlant;

    model JenaSecClarModTakacs
      "Secondary Clarifier Model based on Takacs prepared for a special plant"

      extends WasteWater.Icons.SecClar;
      extends SecClar.Takacs.Interfaces.ratios;
      package SCP = SecClar.Takacs;
      import SI = Modelica.SIunits;
      package WI = Interfaces;
      package WWU = WasteWater.WasteWaterUnits;
      parameter SI.Length hsc=4.0 "height of secondary clarifier";
      parameter Integer n=10 "number of layers of SC model";
      parameter SI.Length zm=hsc/(1.0*n)
        "height of m-th secondary clarifier layer";
      parameter SI.Area Asc=1500.0 "area of secondary clarifier";
      parameter WWU.MassConcentration Xt=3000.0 "threshold for X";
      // total sludge concentration in clarifier feed
      WWU.MassConcentration Xf;
      // layers 1 to 10
      WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},
                {-90,24}})));
      WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,
                47},{112,67}})));
      WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
                -106},{-20,-86}})));
      WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,
                -106},{40,-86}})));
      SCP.bottom_layer S1(
        zm=zm,
        Asc=Asc,
        Xf=Xf,
        rXs=rXs,
        rXbh=rXbh,
        rXba=rXba,
        rXp=rXp,
        rXi=rXi,
        rXnd=rXnd) annotation (Placement(transformation(extent={{-35,-93},{35,-78}})));
      SCP.lower_layer S2(
        zm=zm,
        Asc=Asc,
        Xf=Xf) annotation (Placement(transformation(extent={{-35,-74},{35,-59}})));
      SCP.lower_layer S3(
        zm=zm,
        Asc=Asc,
        Xf=Xf) annotation (Placement(transformation(extent={{-35,-55},{35,-40}})));
      SCP.feed_layer S4(
        zm=zm,
        Asc=Asc,
        Xf=Xf) annotation (Placement(transformation(extent={{-36,-36},{34,-22}})));
      SCP.upper_layer S5(
        zm=zm,
        Asc=Asc,
        Xf=Xf,
        Xt=Xt) annotation (Placement(transformation(extent={{-36,-16},{34,-4}})));
      SCP.upper_layer S6(
        zm=zm,
        Asc=Asc,
        Xf=Xf,
        Xt=Xt) annotation (Placement(transformation(extent={{-36,2},{34,16}})));
      SCP.upper_layer S7(
        zm=zm,
        Asc=Asc,
        Xf=Xf,
        Xt=Xt) annotation (Placement(transformation(extent={{-35,21},{35,36}})));
      SCP.upper_layer S8(
        zm=zm,
        Asc=Asc,
        Xt=Xt,
        Xf=Xf) annotation (Placement(transformation(extent={{-35,40},{35,55}})));
      SCP.upper_layer S9(
        zm=zm,
        Asc=Asc,
        Xf=Xf,
        Xt=Xt) annotation (Placement(transformation(extent={{-35,59},{35,74}})));
      SCP.top_layer S10(
        zm=zm,
        Asc=Asc,
        Xf=Xf,
        Xt=Xt,
        rXs=rXs,
        rXbh=rXbh,
        rXba=rXba,
        rXp=rXp,
        rXi=rXi,
        rXnd=rXnd) annotation (Placement(transformation(extent={{-35,78},{35,93}})));
    equation

      connect(S1.Up, S2.Dn) annotation (Line(points={{-2.22045e-15,-78},{
              -2.22045e-15,-74}}));
      connect(S2.Up, S3.Dn) annotation (Line(points={{-2.22045e-15,-59},{
              -2.22045e-15,-55}}));
      connect(S7.Up, S8.Dn) annotation (Line(points={{-2.22045e-15,36},{
              -2.22045e-15,40}}));
      connect(S9.Up, S10.Dn) annotation (Line(points={{-2.22045e-15,74},{
              -2.22045e-15,78}}));
      connect(S8.Up, S9.Dn) annotation (Line(points={{-2.22045e-15,55},{
              -2.22045e-15,59}}));
      connect(S1.PQw, Waste) annotation (Line(points={{17.5,-93},{30,-93},{30,
              -100}}));
      connect(S10.Out, Effluent) annotation (Line(points={{35,85.5},{67.5,85.5},{
              67.5,57},{100,57}}));
      connect(S1.PQr, Return) annotation (Line(points={{-21,-93},{-30,-93},{-30,
              -100}}));
      connect(S4.Dn, S3.Up) annotation (Line(points={{0,-36},{0,-40}}));
      connect(S4.Up, S5.Dn) annotation (Line(points={{-2,-22},{-2,-16}}));
      connect(S5.Up, S6.Dn) annotation (Line(points={{0,-4},{0,2}}));
      connect(S6.Up, S7.Dn) annotation (Line(points={{0,16},{0,21}}));
      connect(Feed, S4.In) annotation (Line(points={{-100,10},{-67.5,10},{-67.5,
              -28.72},{-35,-28.72}}));

      // total sludge concentration in clarifier feed
      Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

      // ratios of solid components
      rXs = Feed.Xs/Xf;
      rXbh = Feed.Xbh/Xf;
      rXba = Feed.Xba/Xf;
      rXp = Feed.Xp/Xf;
      rXi = Feed.Xi/Xf;
      rXnd = Feed.Xnd/Xf;

      annotation (
        Documentation(info="This component models an ASM1 10 - layer secondary clarifier with 6 layers above the feed_layer (including top_layer)
and 3 layers below the feed_layer (including bottom_layer).

Parameters:
  hsc -  height of clarifier [m]
  n   -  number of layers
  Asc -  surface area of sec. clar. [m2]
  Xt  -  threshold value for Xtss [mg/l]
"));
    end JenaSecClarModTakacs;
    annotation (
      Documentation(info="This package contains example ASM1 wastewater treatment plant models to demonstrate the usage of
the WasteWater.ASM1 library.
Open the models and simulate them according to the description provided in the models.

The following demo models are present:

 - SmallPlant
 - BenchPlant
 - ComplexPlant

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2003, Gerald Reichl
"));
  end Examples;

  package Interfaces
    "Connectors and partial ASM1 model for Wastewater Treatment Modelling"

    extends Modelica.Icons.Library;

    connector WWFlowAsm1in "Inflow connector of ASM1 components"
      package WWU = WasteWater.WasteWaterUnits;
      flow WWU.VolumeFlowRate Q;
      WWU.MassConcentration Si;
      WWU.MassConcentration Ss;
      WWU.MassConcentration Xi;
      WWU.MassConcentration Xs;
      WWU.MassConcentration Xbh;
      WWU.MassConcentration Xba;
      WWU.MassConcentration Xp;
      WWU.MassConcentration So;
      WWU.MassConcentration Sno;
      WWU.MassConcentration Snh;
      WWU.MassConcentration Snd;
      WWU.MassConcentration Xnd;
      WWU.Alkalinity Salk;

      annotation (
        Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics={Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={0,0,191},
              fillPattern=FillPattern.Solid), Text(
              extent={{-88,92},{88,-94}},
              lineColor={255,255,255},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              textString=
                   "%name")}),
        Documentation(info="Connectors WWFlowAsm1in and WWFlowAsm1out are nearly identical.
The difference is in the icons to more easily identify the inflow and outflow
side of a component.
The connector consists of one flow variable and 13 potential variables (ASM1 concentrations).
"));

    end WWFlowAsm1in;

    connector WWFlowAsm1out "Outflow connector of ASM1 components"
      package WWU = WasteWater.WasteWaterUnits;
      flow WWU.VolumeFlowRate Q;
      WWU.MassConcentration Si;
      WWU.MassConcentration Ss;
      WWU.MassConcentration Xi;
      WWU.MassConcentration Xs;
      WWU.MassConcentration Xbh;
      WWU.MassConcentration Xba;
      WWU.MassConcentration Xp;
      WWU.MassConcentration So;
      WWU.MassConcentration Sno;
      WWU.MassConcentration Snh;
      WWU.MassConcentration Snd;
      WWU.MassConcentration Xnd;
      WWU.Alkalinity Salk;
      annotation (
        Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics={Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-88,92},{94,-92}}, textString=
                             "%name")}),
        Documentation(info="Connectors WWFlowAsm1in and WWFlowAsm1out are nearly identical.
The difference is in the icons to more easily identify the inflow and outflow
side of a component.
The connector consists of one flow variable and 13 potential variables (ASM1 concentrations).
"));
    end WWFlowAsm1out;

    connector AirFlow "Airflow connector"
      package WWU = WasteWater.WasteWaterUnits;

      flow WWU.VolumeFlowRate Q_air;
      annotation (
        Documentation(info="The Airflow connector consists of a flow variable describing the exchange of
air between blower and nitrification tank."));
    end AirFlow;

    partial model stoichiometry "ASM1 stoichiometric coefficients"
      // Stoichiometric parameters based on the original ASM1 publication//
      parameter Real Y_h=0.67
        "Heterotrophic Yield [g Xbh COD formed/(g COD utilised)]";
      parameter Real Y_a=0.24
        "Autotrophic Yield [g Xba COD formed/(g N utilised)]";
      parameter Real f_p=0.08 "Fraction of biomass to particulate products [-]";
      parameter Real i_xb=0.08 "Fraction nitrogen in biomass [g N/(g COD)]";
      parameter Real i_xp=0.06
        "Fraction nitrogen in particulate products [g N/(g COD)]";
      annotation (
        Documentation(info=
              "This is a partial model providing the stoichiometric coefficients of the ASM1 model."));
    end stoichiometry;

    partial model ASM1base "Base class of WWTP modelling by ASM1"
      extends Interfaces.stoichiometry;
      package WWU = WasteWater.WasteWaterUnits;

      // parameters based on the original ASM1 publication based on 15 deg C
      Real mu_h "Maximum heterotrophic growth rate f(T) [day^-1]";
      Real b_h "Heterotrophic decay rate f(T) [day^-1]";
      Real mu_a "Maximum autotrophic growth rate f(T) [day^-1]";
      //Real K_nh "Half-saturation (auto. growth) f(T) [g NH-N/m3]";
      Real b_a "Autotrophic decay rate f(T) [day^-1]";
      Real k_a "Ammonification rate f(T) [m3/(g COD day)]";
      Real k_h "Maximum specific hydrolysis rate f(T) [g Xs/(g Xbh COD day)]";
      Real K_x "Half-saturation (hydrolysis) f(T) [g Xs/(g Xbh COD)]";
      parameter Real mu_h_T=4.0
        "Maximum heterotrophic growth rate at T=15 deg C [day^-1]";
      parameter Real b_h_T=0.3
        "Heterotrophic decay rate at T=15 deg C [day^-1]";
      parameter Real mu_a_T=0.5
        "Maximum autotrophic growth rate at T=15 deg C[day^-1]";
      parameter Real b_a_T=0.05 "Autotrophic decay rate at T=15 deg C [day^-1]";
      parameter Real k_a_T=0.05
        "Ammonification rate at T=15 deg C [m3/(g COD day)]";
      parameter Real k_h_T=3.0
        "Maximum specific hydrolysis rate at T=15 deg C [g Xs/(g Xbh COD day)]";
      parameter Real K_x_T=0.1
        "Half-saturation (hydrolysis) at T=15 deg C [g Xs/(g Xbh COD)]";
      parameter Real K_nh=1.0 "Half-saturation (auto. growth) [g NH-N/m3]";
      parameter Real K_s=10.0 "Half-saturation (hetero. growth) [g COD/m3]";
      parameter Real K_oh=0.2 "Half-saturation (hetero. oxygen) [g O/m3]";
      parameter Real K_no=0.5 "Half-saturation (nitrate) [g NO-N/m3]";
      parameter Real K_oa=0.4 "Half-saturation (auto. oxygen) [g O/m3]";
      parameter Real ny_g=0.8 "Anoxic growth rate correction factor [-]";
      parameter Real ny_h=0.8 "Anoxic hydrolysis rate correction factor [-]";
      WWU.MassConcentration Si(fixed=true) "Soluble inert organic matter";
      WWU.MassConcentration Ss(fixed=true) "Readily biodegradable substrate";
      WWU.MassConcentration Xi(fixed=true) "Particulate inert organic matter";
      WWU.MassConcentration Xs(fixed=true) "Slowly biodegradable substrate";
      WWU.MassConcentration Xbh(fixed=true) "Active heterotrophic biomass";
      WWU.MassConcentration Xba(fixed=true) "Active autotrophic biomass";
      WWU.MassConcentration Xp(fixed=true)
        "Particulate products from biomass decay";
      WWU.MassConcentration So(fixed=true) "Dissolved oxygen";
      WWU.MassConcentration Sno(fixed=true) "Nitrate and nitrite nitrogen";
      WWU.MassConcentration Snh(fixed=true) "Ammonium nitrogen";
      WWU.MassConcentration Snd(fixed=true)
        "Soluble biodegradable organic nitrogen";
      WWU.MassConcentration Xnd(fixed=true)
        "Particulate biodegradable organic nitrogen";
      WWU.Alkalinity Salk(fixed=true) "Alkalinity";
      Real p1;
      Real p2;
      Real p3;
      Real p4;
      Real p5;
      Real p6;
      Real p7;
      Real p8;
      Real r1;
      Real r2;
      Real r3;
      Real r4;
      Real r5;
      Real r6;
      Real r7;
      Real r8;
      Real r9;
      Real r10;
      Real r11;
      Real r12;
      Real r13;
      Real inputSi;
      Real inputSs;
      Real inputXi;
      Real inputXs;
      Real inputXbh;
      Real inputXba;
      Real inputXp;
      Real inputSo;
      Real inputSno;
      Real inputSnh;
      Real inputSnd;
      Real inputXnd;
      Real inputSalk;
      Real aeration;

      Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{
                -110,-10},{-90,10}})));
      Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{
                90,-10},{110,10}})));
      Interfaces.WWFlowAsm1out MeasurePort annotation (Placement(transformation(
              extent={{50,40},{60,50}})));
      Modelica.Blocks.Interfaces.RealInput T annotation (Placement(transformation(
              extent={{-110,30},{-90,50}})));
    equation

      // Temperature dependent Kinetic parameters based on 15 deg C //
      // may be adapted to 10 or 20 deg C

      mu_h =mu_h_T;
      b_h =b_h_T;
      mu_a =mu_a_T;
      //K_nh=1.0*exp(0.069*(T.signal[1]-15));
      b_a =b_a_T;
      k_a =k_a_T;
      k_h =k_h_T;
      K_x =K_x_T;

      // Process Rates //

      p1 = mu_h*(Ss/(K_s + Ss))*(So/(K_oh + So))*Xbh;
      p2 = mu_h*(Ss/(K_s + Ss))*(K_oh/(K_oh + So))*(Sno/(K_no + Sno))*ny_g*Xbh;
      p3 = mu_a*(Snh/(K_nh + Snh))*(So/(K_oa + So))*Xba;
      p4 = b_h*Xbh;
      p5 = b_a*Xba;
      p6 = k_a*Snd*Xbh;
      p7 = k_h*((Xs/Xbh)/(K_x + (Xs/Xbh)))*((So/(K_oh + So)) + ny_h*(K_oh/(K_oh
         + So))*(Sno/(K_no + Sno)))*Xbh;
      p8 = p7*Xnd/Xs;

      // biochemical reactions

      r1 = 0;
      r2 = (-p1 - p2)/Y_h + p7;
      r3 = 0;
      r4 = (1 - f_p)*(p4 + p5) - p7;
      r5 = p1 + p2 - p4;
      r6 = p3 - p5;
      r7 = f_p*(p4 + p5);
      r8 = -((1 - Y_h)/Y_h)*p1 - ((4.57 - Y_a)/Y_a)*p3;
      r9 = -((1 - Y_h)/(2.86*Y_h))*p2 + p3/Y_a;
      r10 = -i_xb*(p1 + p2) - (i_xb + (1/Y_a))*p3 + p6;
      r11 = -p6 + p8;
      r12 = (i_xb - f_p*i_xp)*(p4 + p5) - p8;
      r13 = -i_xb/14*p1 + ((1 - Y_h)/(14*2.86*Y_h) - (i_xb/14))*p2 - ((i_xb/14)
         + 1/(7*Y_a))*p3 + p6/14;

      // derivatives

      der(Si) = inputSi + r1;
      der(Ss) = inputSs + r2;
      der(Xi) = inputXi + r3;
      der(Xs) = inputXs + r4;
      der(Xbh) = inputXbh + r5;
      der(Xba) = inputXba + r6;
      der(Xp) = inputXp + r7;
      der(So) = inputSo + r8 + aeration;
      der(Sno) = inputSno + r9;
      der(Snh) = inputSnh + r10;
      der(Snd) = inputSnd + r11;
      der(Xnd) = inputXnd + r12;
      der(Salk) = inputSalk + r13;

      // Outputs

      Out.Q + In.Q = 0;
      Out.Si = Si;
      Out.Ss = Ss;
      Out.Xi = Xi;
      Out.Xs = Xs;
      Out.Xbh = Xbh;
      Out.Xba = Xba;
      Out.Xp = Xp;
      Out.So = So;
      Out.Sno = Sno;
      Out.Snh = Snh;
      Out.Snd = Snd;
      Out.Xnd = Xnd;
      Out.Salk = Salk;

      MeasurePort.Si = Si;
      MeasurePort.Ss = Ss;
      MeasurePort.Xi = Xi;
      MeasurePort.Xs = Xs;
      MeasurePort.Xbh = Xbh;
      MeasurePort.Xba = Xba;
      MeasurePort.Xp = Xp;
      MeasurePort.So = So;
      MeasurePort.Sno = Sno;
      MeasurePort.Snh = Snh;
      MeasurePort.Snd = Snd;
      MeasurePort.Xnd = Xnd;
      MeasurePort.Salk = Salk;

      annotation (
        Documentation(info="This partial model provides connectors and equations that are needed in the biological
components (nitrification and denitrification tank) for ASM1 wastewater treatment plant models.
Parameters are coded according the ASM1 [1,2] standard distribution.
Changes to this parameters are subject to the modeller.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de


References:

[1] M. Henze and C.P.L. Grady Jr and W. Gujer and G.v.R. Marais and T. Matsuo:
    Activated Sludge Model No.1. IAWQ, 1987.
[2] M. Henze and W.Gujer and T. Mino and. M.v. Loosdrecht: Activated Sludge
    Models ASM1, ASM2, ASM2d, and ASM3. IWA Task Group on Mathematical Modelling
    for Design and Operation of Biological Wastewater Treatment, 2000.


This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2002, Gerald Reichl
"));
    end ASM1base;
    annotation (
      Documentation(info="This package contains connectors and interfaces (partial models) for
wastewater treatment components based on the Acticated Sludge Model No.1 (ASM1).

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2001, Gerald Reichl
"));
  end Interfaces;

  package PreClar "Primary clarifier modelling based on ASM1"

    extends Modelica.Icons.Library;

    model preclar1 "Dynamic ASM1 Primary Clarifier Model"
      import Modelica.Math.log;

      // dynamic primary clarifier tank, based on Otterpohl
      // to be used for feed forward calculation, e.g. influent data needed

      package WWU = WasteWaterUnits;
      extends WasteWater.Icons.preclar1;

      // tank specific parameters

      parameter Modelica.SIunits.Volume V=500
        "Volume of primary clarifier tank";
      Real hrt_h "hydraulic residence time in primary sedimentation tank [h]";

        // Real hrt_min "hydraulic residence time in primary sedimentation tank [min]";
      Real n_COD "efficiency of COD removal [%]";
      Real n_X "efficiency transformed to particulate fractions [%]";
      WWU.MassConcentration Si "Soluble inert organic matter";
      WWU.MassConcentration Ss "Readily biodegradable substrate";
      WWU.MassConcentration Xi "Particulate inert organic matter";
      WWU.MassConcentration Xs "Slowly biodegradable substrate";
      WWU.MassConcentration Xbh "Active heterotrophic biomass";
      WWU.MassConcentration Xba "Active autotrophic biomass";
      WWU.MassConcentration Xp "Particulate products from biomass decay";
      WWU.MassConcentration So "Dissolved oxygen";
      WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
      WWU.MassConcentration Snh "Ammonium nitrogen";
      WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
      WWU.MassConcentration Xnd "Particulate biodegradable organic nitrogen";
      WWU.Alkalinity Salk "Alkalinity";
      Real CODin;
      Real CODout;
      Real XCODin;
      Real H;
      BSM1.Interfaces.WWFlowAsm1in In
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
      BSM1.Interfaces.WWFlowAsm1out Out
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      BSM1.Interfaces.WWFlowAsm1out MeasurePort
        annotation (Placement(transformation(extent={{32,90},{42,100}})));
    equation

      // calculation of the hydraulic residence time
      hrt_h = V/In.Q*24;
      //hrt_min = V/In.Q * 24 * 60;

        // n_COD according Otterpohl and Freund 1992 "Dynamic Models for Clarifiers"
      n_COD = 2.7*(log(hrt_h*hrt_h) + 9)/100;
      // n_COD according Otterpohl 1995, Dissertation
      // n_COD = (1.45 + 6.15 * log(hrt_min))/100;

      XCODin = In.Xi + In.Xs + In.Xbh + In.Xba + In.Xp;
      // particulate COD in the influent
      CODin = In.Si + In.Ss + XCODin;
      // total COD in the influent

      CODout = Out.Si + Out.Ss + Out.Xi + Out.Xs + Out.Xbh + Out.Xba + Out.Xp;

      H = n_COD*CODin/XCODin;

      // n_X can not be greater than 1
      // therefore is this check
      n_X = if H > 0.95 then 0.95 else if H < 0.05 then 0.05 else H;

      // in this case the model needs to be modified by a new n_COD
      // n_COD_? = (2.88*XCODin/CODin - 0.118) * n_COD;

      // volume dependent dilution term of each concentration

      der(Si) = (In.Si - Si)*In.Q/V;
      der(Ss) = (In.Ss - Ss)*In.Q/V;
      der(Xi) = (In.Xi - Xi)*In.Q/V;
      der(Xs) = (In.Xs - Xs)*In.Q/V;
      der(Xbh) = (In.Xbh - Xbh)*In.Q/V;
      der(Xba) = (In.Xba - Xba)*In.Q/V;
      der(Xp) = (In.Xp - Xp)*In.Q/V;
      der(So) = (In.So - So)*In.Q/V;
      der(Sno) = (In.Sno - Sno)*In.Q/V;
      der(Snh) = (In.Snh - Snh)*In.Q/V;
      der(Snd) = (In.Snd - Snd)*In.Q/V;
      der(Xnd) = (In.Xnd - Xnd)*In.Q/V;
      der(Salk) = (In.Salk - Salk)*In.Q/V;

      // Outputs
      // this is just a reduction of particulate substances; n_X*X is not stored
      // so the amount of primary sludge removed is not calculated
      Out.Q + In.Q = 0;

      Out.Si = Si;
      Out.Ss = Ss;
      Out.Xi = (1 - n_X)*Xi;
      Out.Xs = (1 - n_X)*Xs;
      Out.Xbh = (1 - n_X)*Xbh;
      Out.Xba = (1 - n_X)*Xba;
      Out.Xp = (1 - n_X)*Xp;
      Out.So = So;
      Out.Sno = Sno;
      Out.Snh = Snh;
      Out.Snd = Snd;
      Out.Xnd = (1 - n_X)*Xnd;
      Out.Salk = Salk;

      MeasurePort.Si = Si;
      MeasurePort.Ss = Ss;
      MeasurePort.Xi = (1 - n_X)*Xi;
      MeasurePort.Xs = (1 - n_X)*Xs;
      MeasurePort.Xbh = (1 - n_X)*Xbh;
      MeasurePort.Xba = (1 - n_X)*Xba;
      MeasurePort.Xp = (1 - n_X)*Xp;
      MeasurePort.So = So;
      MeasurePort.Sno = Sno;
      MeasurePort.Snh = Snh;
      MeasurePort.Snd = Snd;
      MeasurePort.Xnd = (1 - n_X)*Xnd;
      MeasurePort.Salk = Salk;
      annotation (
        Documentation(info="This is an ASM1 dynamic primary clarifier model based on the theory
by Otterpohl and Freund.

Parameter:
  V - volume of the preclarifier
"));
    end preclar1;

    model preclar2 "Static ASM1 Primary Clarifier Model"
      // static primary clarifier tank, based on Otterpohl
      // to be used for feed forward calculation, e.g. influent data needed

      import Modelica.Math.log;

      package WWU = WasteWaterUnits;
      extends WasteWater.Icons.preclar2;

      // tank specific parameters

      parameter Modelica.SIunits.Volume V=500
        "Volume of primary clarifier tank";
      Real hrt_h "hydraulic residence time in primary sedimentation tank [h]";

        //Real hrt_min "hydraulic residence time in primary sedimentation tank [min]";
      Real n_COD "efficiency of COD removal [%]";
      Real n_X "efficiency transformed to particulate fractions [%]";

      WWU.MassConcentration Si "Soluble inert organic matter";
      WWU.MassConcentration Ss "Readily biodegradable substrate";
      WWU.MassConcentration Xi "Particulate inert organic matter";
      WWU.MassConcentration Xs "Slowly biodegradable substrate";
      WWU.MassConcentration Xbh "Active heterotrophic biomass";
      WWU.MassConcentration Xba "Active autotrophic biomass";
      WWU.MassConcentration Xp "Particulate products from biomass decay";
      WWU.MassConcentration So "Dissolved oxygen";
      WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
      WWU.MassConcentration Snh "Ammonium nitrogen";
      WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
      WWU.MassConcentration Xnd "Particulate biodegradable organic nitrogen";
      WWU.Alkalinity Salk "Alkalinity";

      Real CODin;
      Real CODout;
      Real XCODin;
      Real H;
      BSM1.Interfaces.WWFlowAsm1in In
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
      BSM1.Interfaces.WWFlowAsm1out Out
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      BSM1.Interfaces.WWFlowAsm1out MeasurePort
        annotation (Placement(transformation(extent={{32,90},{42,100}})));
    equation

      // calculation of the hydraulic residence time
      hrt_h = V/In.Q*24;
      //hrt_min = V/In.Q * 24 * 60;

        // n_COD according Otterpohl and Freund 1992 "Dynamic Models for Clarifiers"
      n_COD = 2.7*(log(hrt_h*hrt_h) + 9)/100;
      // n_COD according Otterpohl 1995, Dissertation
      // n_COD = (1.45 + 6.15 * log(hrt_min))/100;

      XCODin = In.Xi + In.Xs + In.Xbh + In.Xba + In.Xp;
      // particulate COD in the influent
      CODin = In.Si + In.Ss + XCODin;
      // total COD in the influent

      CODout = Out.Si + Out.Ss + Out.Xi + Out.Xs + Out.Xbh + Out.Xba + Out.Xp;

      H = n_COD*CODin/XCODin;
      // n_X can not be greater than 1
      // therefore is this check
      n_X = if H > 0.95 then 0.95 else if H < 0.05 then 0.05 else H;
      // in this case the model needs to be modified by a new n_COD
      // n_COD_? = (2.88*XCODin/CODin - 0.118) * n_COD;

      // volume dependent dilution term of each concentration

      0 = (In.Si - Si)*In.Q/V;
      0 = (In.Ss - Ss)*In.Q/V;
      0 = (In.Xi - Xi)*In.Q/V;
      0 = (In.Xs - Xs)*In.Q/V;
      0 = (In.Xbh - Xbh)*In.Q/V;
      0 = (In.Xba - Xba)*In.Q/V;
      0 = (In.Xp - Xp)*In.Q/V;
      0 = (In.So - So)*In.Q/V;
      0 = (In.Sno - Sno)*In.Q/V;
      0 = (In.Snh - Snh)*In.Q/V;
      0 = (In.Snd - Snd)*In.Q/V;
      0 = (In.Xnd - Xnd)*In.Q/V;
      0 = (In.Salk - Salk)*In.Q/V;

      // Outputs
      // this is just a reduction of particulate substances; n_X*X is not stored
      // so the amount of primary sludge removed is not calculated
      Out.Q + In.Q = 0;
      Out.Si = Si;
      Out.Ss = Ss;
      Out.Xi = (1 - n_X)*Xi;
      Out.Xs = (1 - n_X)*Xs;
      Out.Xbh = (1 - n_X)*Xbh;
      Out.Xba = (1 - n_X)*Xba;
      Out.Xp = (1 - n_X)*Xp;
      Out.So = So;
      Out.Sno = Sno;
      Out.Snh = Snh;
      Out.Snd = Snd;
      Out.Xnd = (1 - n_X)*Xnd;
      Out.Salk = Salk;

      MeasurePort.Si = Si;
      MeasurePort.Ss = Ss;
      MeasurePort.Xi = (1 - n_X)*Xi;
      MeasurePort.Xs = (1 - n_X)*Xs;
      MeasurePort.Xbh = (1 - n_X)*Xbh;
      MeasurePort.Xba = (1 - n_X)*Xba;
      MeasurePort.Xp = (1 - n_X)*Xp;
      MeasurePort.So = So;
      MeasurePort.Sno = Sno;
      MeasurePort.Snh = Snh;
      MeasurePort.Snd = Snd;
      MeasurePort.Xnd = (1 - n_X)*Xnd;
      MeasurePort.Salk = Salk;

      annotation (
        Documentation(info="This is an ASM1 static primary clarifier model based on the theory
by Otterpohl and Freund.

Parameter:
  V - volume of the preclarifier

"));
    end preclar2;

    model preclar3 "Inverse ASM1 Static Primary Clarifier Model"
      // static primary clarifier tank

        // to be used for backward calculation, e.g. effluent concentration data needed
      // signals need to be in the secuence COD, Sno, Snh, pH in the inputtable

      import Modelica.Math.log;
      extends WasteWater.Icons.preclar2;

      package WWU = WasteWater.WasteWaterUnits;

        // Interfaces.MeasurePort MeasurePort annotation (extent=[32, 90; 42, 100]);

      // tank specific parameters
      parameter Modelica.SIunits.Volume V=500
        "Volume of primary clarifier tank";
      parameter Real aSi=5/100
        "Fraction of Si of the total COD in the influent";
      parameter Real aSs=15/100
        "Fraction of Ss of the total COD in the influent";
      parameter Real aXi=15/100
        "Fraction of Xi of the total COD in the influent";
      parameter Real aXs=45/100
        "Fraction of Xs of the total COD in the influent";
      parameter Real aXbh=20/100
        "Fraction of Xbh of the total COD in the influent";
      parameter Real aXba=0/100
        "Fraction of Xba of the total COD in the influent";
      parameter Real aXp=0/100
        "Fraction of Xp of the total COD in the influent";
      parameter Real aSo=0.0 "Dissolved oxygen in the inflow [mg/l]";
      parameter Real aSnd=1/100 "Fraction Snd of Ss in the influent";
      parameter Real aXnd=3/100 "Fraction Xnd of Xs in the influent";
      parameter Real n_corr=1.0 "Correction faktor for the efficiency function";

      Real hrt_h "hydraulic residence time in primary sedimentation tank [h]";

        //Real hrt_min "hydraulic residence time in primary sedimentation tank [min]";
      Real n_COD "efficiency of COD removal [%]";
      Real n_X "efficiency transformed to particulate fractions [%]";

      WWU.MassConcentration Si "Soluble inert organic matter";
      WWU.MassConcentration Ss "Readily biodegradable substrate";
      WWU.MassConcentration Xi "Particulate inert organic matter";
      WWU.MassConcentration Xs "Slowly biodegradable substrate";
      WWU.MassConcentration Xbh "Active heterotrophic biomass";
      WWU.MassConcentration Xba "Active autotrophic biomass";
      WWU.MassConcentration Xp "Particulate products from biomass decay";
      WWU.MassConcentration So "Dissolved oxygen";
      WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
      WWU.MassConcentration Snh "Ammonium nitrogen";
      WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
      WWU.MassConcentration Xnd "Particulate biodegradable organic nitrogen";
      WWU.Alkalinity Salk "Alkalinity";
      Real COD;
      Real CODin;
      Real CODout;
      Real XCOD;
      Real H;
      BSM1.Interfaces.WWFlowAsm1in In
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
      BSM1.Interfaces.WWFlowAsm1out Out
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      Modelica.Blocks.Interfaces.RealInput MeasurePort[4]
        annotation (Placement(transformation(
            origin={38,90},
            extent={{-10,-10},{10,10}},
            rotation=270)));
    equation

      // calculation of the hydraulic residence time
      hrt_h = V/In.Q*24;
      //hrt_min = V/In.Q * 24 * 60;

        // n_COD according Otterpohl and Freund 1992 "Dynamic Models for Clarifiers"
      n_COD = n_corr*2.7*(log(hrt_h*hrt_h) + 9)/100;
      // n_COD according Otterpohl 1995, Dissertation
      // n_COD = (1.45 + 6.15 * log(hrt_min))/100;

      XCOD = In.Xi + In.Xs + In.Xbh + In.Xba + In.Xp;
      // particulate COD in the influent
      COD = In.Si + In.Ss + XCOD;
      // total COD in the influent

      CODin =MeasurePort[1]/(1 - n_COD);
      // total COD in the influent
      // above two CODs sould be the same

      CODout = Out.Si + Out.Ss + Out.Xi + Out.Xs + Out.Xbh + Out.Xba + Out.Xp;
      // this should be the same as MeasurePort.signal[1]

      H = n_COD*COD/XCOD;
      // n_X can not be greater than 1
      // therefor this check
      n_X = if H > 0.95 then 0.95 else if H < 0.05 then 0.05 else H;
      // in this case the model needs to be modified by a new n_COD
      // n_COD_? = (2.88*XCODin/CODin - 0.118) * n_COD;

      // volume dependent dilution term of each concentration

      0 = (In.Si - Si)*In.Q/V;
      0 = (In.Ss - Ss)*In.Q/V;
      0 = (In.Xi - Xi)*In.Q/V;
      0 = (In.Xs - Xs)*In.Q/V;
      0 = (In.Xbh - Xbh)*In.Q/V;
      0 = (In.Xba - Xba)*In.Q/V;
      0 = (In.Xp - Xp)*In.Q/V;
      0 = (In.So - So)*In.Q/V;
      0 = (In.Sno - Sno)*In.Q/V;
      0 = (In.Snh - Snh)*In.Q/V;
      0 = (In.Snd - Snd)*In.Q/V;
      0 = (In.Xnd - Xnd)*In.Q/V;
      0 = (In.Salk - Salk)*In.Q/V;

      Out.Q + In.Q = 0;

      // Inputs
      In.Si = aSi*CODin;
      In.Ss = aSs*CODin;
      In.Xi = aXi*CODin;
      In.Xs = aXs*CODin;
      In.Xbh = aXbh*CODin;
      In.Xba = aXba*CODin;
      In.Xp = aXp*CODin;
      In.So = aSo;
      In.Sno =MeasurePort[2];
      In.Snh =MeasurePort[3];
      In.Snd = aSnd*In.Ss;
      In.Xnd = aXnd*In.Xs;
      In.Salk =1.8*exp(MeasurePort[4] - 6.4);

      // Outputs
      Out.Si = Si;
      Out.Ss = Ss;
      Out.Xi = (1 - n_X)*Xi;
      Out.Xs = (1 - n_X)*Xs;
      Out.Xbh = (1 - n_X)*Xbh;
      Out.Xba = (1 - n_X)*Xba;
      Out.Xp = (1 - n_X)*Xp;
      Out.So = So;
      Out.Sno = Sno;
      Out.Snh = Snh;
      Out.Snd = Snd;
      Out.Xnd = (1 - n_X)*Xnd;
      Out.Salk = Salk;

      annotation (
        Documentation(info="This is a special case of the ASM1 static primary clarifier model.
Here measurement data at the end (effluent) of the preclarifier needs to be provided.
This is typical for some real plants. Influent is then calculated.

Parameters:
  V   - volume of the preclarifier
  aS* - fractions of e.g. COD in influent (soluble components)
  aX* - fractions of e.g. COD in influent (particular components)
  n_corr- correction factor for the efficiency function

Dimension of InPort is 4.
  1 - Chemical Oxygen Demand (COD) at effluent of primary clarifier
  2 - nitrate nitrogen (Sno) at effluent of primary clarifier
  3 - ammonium nitrogen (Snh) at effluent of primary clarifier
  4 - pH-value at effluent of primary clarifier

"));
    end preclar3;
    annotation (
      Documentation(info="This package provides one dynamic and two static ASM1 primary clarifier
models based on Otterpohl [1].

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de


Reference:

[1] R. Otterpohl and M. Freund: Dynamic models for clarifier of activated sludge
    plants with dry and wet weather flows. Water Science and Technology. 26 (1992), pp. 1391-1400.

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2001, Gerald Reichl
"));
  end PreClar;

  package SecClar "Library for secondary settling tank modelling based on ASM1"
  extends Modelica.Icons.Library;

    package Haertel "Secondary settling tank modelling by Haertel (ASM1)"
      extends Modelica.Icons.Library;

      package Interfaces
        "Connectors and partial models for ASM1 Secondary Clarifier Model by Haertel"

        extends Modelica.Icons.Library;

        connector UpperLayerPin "Connector above influent layer"

          package WWU = WasteWater.WasteWaterUnits;

          // effluent flow
          flow WWU.VolumeFlowRate Qe;

          // sedimentation flux
          flow WWU.SedimentationFlux SedFlux;

          // total sludge concentration (m-1)-th layer (dn=down)
          WWU.MassConcentration X_dn;

          // soluble components
          WWU.MassConcentration Si;
          WWU.MassConcentration Ss;
          WWU.MassConcentration So;
          WWU.MassConcentration Sno;
          WWU.MassConcentration Snh;
          WWU.MassConcentration Snd;
          WWU.Alkalinity Salk;
          annotation (
            Documentation(info=
                  "Connector for ASM1 information and mass exchange between layers above the influent layer (feed_layer)."));
        end UpperLayerPin;

        connector LowerLayerPin "Connector below influent layer"
          package WWU = WasteWater.WasteWaterUnits;

          // return and waste sludge flow Qr, Qw
          flow WWU.VolumeFlowRate Qr;
          flow WWU.VolumeFlowRate Qw;
          // sedimentation flux
          flow WWU.SedimentationFlux SedFlux;
          // total sludge concentration in m-th layer
          WWU.MassConcentration X;

            // total sludge concentration and sink velocity in (m-1)-th layer (dn=down)
          WWU.MassConcentration X_dn;
          WWU.SedimentationVelocity vS_dn;
          // soluble components
          WWU.MassConcentration Si;
          WWU.MassConcentration Ss;
          WWU.MassConcentration So;
          WWU.MassConcentration Sno;
          WWU.MassConcentration Snh;
          WWU.MassConcentration Snd;
          WWU.Alkalinity Salk;
          annotation (
            Documentation(info=
                  "Connector for ASM1 information and mass exchange between layers below the influent layer (feed_layer)."));
        end LowerLayerPin;

        partial model SCParam "partial model providing clarifier parameters"
          import SI = Modelica.SIunits;
          package WWU = WasteWater.WasteWaterUnits;
          parameter SI.Length zm;
          parameter SI.Area Asc;
          parameter WWU.SludgeVolumeIndex ISV;
          annotation (
            Documentation(info="partial model providing clarifier parameters"));
        end SCParam;

        partial model SCVar "partial models providing variables"
          package WWU = WasteWater.WasteWaterUnits;
          WWU.MassConcentration X "total sludge concentration in m-th layer";
          WWU.SedimentationVelocity vS "sink velocity in m-th layer";
          WWU.SedimentationFlux Jsm "sedimentation flux m-th layer";

          WWU.MassConcentration Si "Soluble inert organic matter";
          WWU.MassConcentration Ss "Readily biodegradable substrate";
          WWU.MassConcentration So "Dissolved oxygen";
          WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
          WWU.MassConcentration Snh "Ammonium nitrogen";
          WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
          WWU.Alkalinity Salk "Alkalinity";
          annotation (
            Documentation(info="partial models providing variables"));
        end SCVar;

        partial model ratios "partial model for ratios of solid components"
          // ratios of solid components
          Real rXi;
          Real rXs;
          Real rXbh;
          Real rXba;
          Real rXp;
          Real rXnd;
          annotation (
            Documentation(info="partial model for ASM1 ratios of solid components"));
        end ratios;

        function vSfun "Sedimentation velocity function"

          // total sludge concentration in m-th layer in g/m3 or mg/l
          input Real X;
          //Sludge Volume Index
          input Real ISV;

          // sink velocity in m/d
          output Real vS;
        protected
          Real v0 "maximum settling velocity";
          Real nv "exponent as part of the Vesilind equation";
        algorithm
          v0 := (17.4*(exp(-0.0113*ISV)) + 3.931)*24;
          //[m/d]
          nv := (-0.9834*(exp(-0.00581*ISV)) + 1.043);
          //[l/g]
          vS := v0*exp(-nv*X/1000);
          annotation (
            Documentation(info="Sedimentation velocity function"));
        end vSfun;

        function omega "Omega correction function by Haertel"
          //vertical coordinate, bottom: z=0
          input Real z;
          // total sludge concentration in clarifier feed
          input Real Xf;
          //height of secondary clarifier
          input Real hsc;
          //height of m-th secondary clarifier layer
          input Real zm;
          //Sludge Volume Index
          input Real ISV;
          //number of layers above feed layer
          input Integer i;
          // correction function omega by Haertel based on [g/l]
          output Real omega;

        protected
          Real Xc "solids concentration at compression point";
          Real nv "exponent as part of the Vesilind equation";
          Real ht "height of transition point";
          Real hc "height of compressing point";
          Real B3;
          Real B4;

        algorithm
          Xc := 480/ISV;
          nv := 1.043 - 0.9834*exp(-0.00581*ISV);
          hc := (Xf/1000)*(hsc - zm*(i + 0.5))/Xc*(1.0 - 1.0/(Xc*nv));
          // unit change
          ht := min(2.0*hc, hsc - zm*(i + 0.5));

          B4 := 1.0 + 2.0*ISV/(100.0 + ISV);
          B3 := -((2*ISV + 100.0)/ISV)*hc^B4;

          omega := (1.0 - B3*ht^(-B4))/(1.0 - B3*z^(-B4));
          omega := min(1.0, omega);
          annotation (
            Documentation(info=
                  "This is Haertels omega correction function for the settling process."));
        end omega;
        annotation (
          Documentation(info="This package contains connectors and interfaces (partial models) for
the ASM1 secondary clarifier model based on Haertel [1] (settling velocity and omega correction function).

References:

[1] L. Haertel: Modellansaetze zur dynamischen Simulation des Belebtschlammverfahrens.
    TH Darmstadt, Dissertation, 1990.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2002, Gerald Reichl
"));
      end Interfaces;

      model SecClarModHaertel
        "ASM1 Secondary Settling Tank Model based on Haertel"
        import WasteWater;

        extends WasteWater.Icons.SecClar;
        extends Interfaces.ratios;
        package SCP = Haertel;
        import SI = Modelica.SIunits;
        package WI = WasteWater.BSM1.Interfaces;
        package WWU = WasteWater.WasteWaterUnits;
        parameter SI.Length hsc=4.0 "height of secondary clarifier";
        parameter Integer n=10 "number of layers of SC model";
        parameter SI.Length zm=hsc/(1.0*n)
          "height of m-th secondary clarifier layer";
        parameter SI.Area Asc=1500.0 "area of secondary clarifier";
        parameter WWU.SludgeVolumeIndex ISV=130 "Sludge Volume Index";
        parameter Integer i=2
          "number of layers above current feed layer in this model";

        // total sludge concentration in clarifier feed
        WWU.MassConcentration Xf;

        // layers 1 to 10
        SCP.bottom_layer S1(
          zm=zm,
          Asc=Asc,
          ISV=ISV,
          rXs=rXs,
          rXbh=rXbh,
          rXba=rXba,
          rXp=rXp,
          rXi=rXi,
          rXnd=rXnd) annotation (Placement(transformation(extent={{-35,-93},{35,-78}})));
        SCP.lower_layer S2(
          hsc=hsc,
          zm=zm,
          z=(zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-74},{35,-59}})));
        SCP.lower_layer S3(
          hsc=hsc,
          zm=zm,
          z=(2*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-55},{35,-40}})));
        SCP.lower_layer S4(
          hsc=hsc,
          zm=zm,
          z=(3*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-36},{35,-21}})));
        SCP.lower_layer S5(
          hsc=hsc,
          zm=zm,
          z=(4*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-17},{35,-2}})));
        SCP.lower_layer S6(
          hsc=hsc,
          zm=zm,
          z=(5*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,2},{35,17}})));
        SCP.lower_layer S7(
          hsc=hsc,
          zm=zm,
          z=(6*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,21},{35,36}})));
        SCP.feed_layer S8(
          hsc=hsc,
          zm=zm,
          z=(7*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,40},{35,55}})));
        SCP.upper_layer S9(
          zm=zm,
          Asc=Asc,
          ISV=ISV) annotation (Placement(transformation(extent={{-35,59},{35,74}})));
        SCP.top_layer S10(
          zm=zm,
          Asc=Asc,
          ISV=ISV,
          rXs=rXs,
          rXbh=rXbh,
          rXba=rXba,
          rXp=rXp,
          rXi=rXi,
          rXnd=rXnd) annotation (Placement(transformation(extent={{-35,78},{35,93}})));
        WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},
                  {-90,24}})));
        WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,
                  47},{112,67}})));
        WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
                  -106},{-20,-86}})));
        WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,
                  -106},{40,-86}})));
      equation

        connect(S1.Up, S2.Dn) annotation (Line(points={{-2.22045e-15,-78},{
                -2.22045e-15,-74}}));
        connect(S2.Up, S3.Dn) annotation (Line(points={{-2.22045e-15,-59},{
                -2.22045e-15,-55}}));
        connect(S3.Up, S4.Dn) annotation (Line(points={{-2.22045e-15,-40},{
                -2.22045e-15,-36}}));
        connect(S5.Up, S6.Dn) annotation (Line(points={{-2.22045e-15,-2},{
                -2.22045e-15,2}}));
        connect(S6.Up, S7.Dn) annotation (Line(points={{-2.22045e-15,17},{
                -2.22045e-15,21}}));
        connect(S7.Up, S8.Dn) annotation (Line(points={{-2.22045e-15,36},{
                -2.22045e-15,40}}));
        connect(S9.Up, S10.Dn) annotation (Line(points={{-2.22045e-15,74},{
                -2.22045e-15,78}}));
        connect(S4.Up, S5.Dn) annotation (Line(points={{-2.22045e-15,-21},{
                -2.22045e-15,-17}}));
        connect(S8.Up, S9.Dn) annotation (Line(points={{-2.22045e-15,55},{
                -2.22045e-15,59}}));
        connect(Feed, S8.In) annotation (Line(points={{-98,14},{-98,47.8},{-33,47.8}}));
        connect(S1.PQw, Waste) annotation (Line(points={{17.5,-93},{17.5,-100},{30,
                -100}}));
        connect(S10.Out, Effluent) annotation (Line(points={{35,85.5},{67.5,85.5},{
                67.5,57},{100,57}}));
        connect(S1.PQr, Return) annotation (Line(points={{-21,-93},{-21,-100},{-30,
                -100}}));

        // total sludge concentration in clarifier feed
        Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

        // ratios of solid components
        rXs = Feed.Xs/Xf;
        rXbh = Feed.Xbh/Xf;
        rXba = Feed.Xba/Xf;
        rXp = Feed.Xp/Xf;
        rXi = Feed.Xi/Xf;
        rXnd = Feed.Xnd/Xf;

        annotation (
          Documentation(info="This component models an ASM1 10 - layer secondary clarifier with 4 layers above the feed_layer (including top_layer)
and 5 layers below the feed_layer (including bottom_layer) based on Haertel`s theory.

Parameters:
  hsc -  height of clarifier [m]
  n   -  number of layers
  Asc -  surface area of sec. clar. [m2]
  ISV -  Sludge Volume Index [ml/g]
  i   -  number of layers above feed layer
"));
      end SecClarModHaertel;

      model bottom_layer "Bottom layer of Haertel`s SC model"
        import WWSC = WasteWater.BSM1.SecClar.Haertel.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        extends WWSC.ratios;

        BSM1.Interfaces.WWFlowAsm1out PQr
          annotation (Placement(transformation(extent={{-70,-110},{-50,-90}})));
        BSM1.Interfaces.WWFlowAsm1out PQw
          annotation (Placement(transformation(extent={{40,-110},{60,-90}})));
        WWSC.LowerLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, ISV);

        // sedimentation flux in bottom layer
        Jsm = 0.0;

        // ODE of solid component
        der(X) = ((Up.Qr + Up.Qw)/Asc*(Up.X - X) + Up.SedFlux)/zm;

        // ODEs of soluble components
        der(Si) = (Up.Qr + Up.Qw)*(Up.Si - Si)/(Asc*zm);
        der(Ss) = (Up.Qr + Up.Qw)*(Up.Ss - Ss)/(Asc*zm);
        der(So) = (Up.Qr + Up.Qw)*(Up.So - So)/(Asc*zm);
        der(Sno) = (Up.Qr + Up.Qw)*(Up.Sno - Sno)/(Asc*zm);
        der(Snh) = (Up.Qr + Up.Qw)*(Up.Snh - Snh)/(Asc*zm);
        der(Snd) = (Up.Qr + Up.Qw)*(Up.Snd - Snd)/(Asc*zm);
        der(Salk) = (Up.Qr + Up.Qw)*(Up.Salk - Salk)/(Asc*zm);

        // upward connection
        Up.vS_dn = vS;
        Up.X_dn = X;

        // return and waste sludge volume flow rates
        PQr.Q + Up.Qr = 0;
        PQw.Q + Up.Qw = 0;

        // return sludge flow, solid and soluble components (ASM1)
        PQr.Si = Si;
        PQr.Ss = Ss;
        PQr.Xi = rXi*X;
        PQr.Xs = rXs*X;
        PQr.Xbh = rXbh*X;
        PQr.Xba = rXba*X;
        PQr.Xp = rXp*X;
        PQr.So = So;
        PQr.Sno = Sno;
        PQr.Snh = Snh;
        PQr.Snd = Snd;
        PQr.Xnd = rXnd*X;
        PQr.Salk = Salk;

        // waste sludge flow, solid and soluble components (ASM1)
        PQw.Si = Si;
        PQw.Ss = Ss;
        PQw.Xi = rXi*X;
        PQw.Xs = rXs*X;
        PQw.Xbh = rXbh*X;
        PQw.Xba = rXba*X;
        PQw.Xp = rXp*X;
        PQw.So = So;
        PQw.Sno = Sno;
        PQw.Snh = Snh;
        PQw.Snd = Snd;
        PQw.Xnd = rXnd*X;
        PQw.Salk = Salk;

        annotation (
          Documentation(info="This class models the lowest layer of an ASM1 secondary clarifier based on Haertel.

No sedimentation flux (mass exchange) with underneath but hydraulic and sedimentation flux (same direction)
with above layer.
From here return and waste sludge is removed.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end bottom_layer;

      model lower_layer "Layer below influent of Haertel`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Haertel.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWU.MassConcentration Xf "sludge concentration in clarifier feed";
        SI.Length z "vertical coordinate of current layer";

        parameter SI.Length hsc;
        parameter Integer i "number of layers above feed layer";
        Real omega;

        WWSC.LowerLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
        WWSC.LowerLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, ISV);
        omega = WWSC.omega(z, Xf, hsc, zm, ISV, i);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = if vS < Dn.vS_dn then omega*(vS*X) else omega*min(vS*X, Dn.vS_dn*Dn.
          X_dn);

        // ODE of solid component
        der(X) = ((Up.Qr + Up.Qw)/Asc*(Up.X - X) + Up.SedFlux - Jsm)/zm;

        // ODEs of soluble components
        der(Si) = (Up.Qr + Up.Qw)*(Up.Si - Si)/(Asc*zm);
        der(Ss) = (Up.Qr + Up.Qw)*(Up.Ss - Ss)/(Asc*zm);
        der(So) = (Up.Qr + Up.Qw)*(Up.So - So)/(Asc*zm);
        der(Sno) = (Up.Qr + Up.Qw)*(Up.Sno - Sno)/(Asc*zm);
        der(Snh) = (Up.Qr + Up.Qw)*(Up.Snh - Snh)/(Asc*zm);
        der(Snd) = (Up.Qr + Up.Qw)*(Up.Snd - Snd)/(Asc*zm);
        der(Salk) = (Up.Qr + Up.Qw)*(Up.Salk - Salk)/(Asc*zm);

        // downward connections
        Dn.Qr + Up.Qr = 0;
        Dn.Qw + Up.Qw = 0;

        Dn.X = X;
        Dn.SedFlux = -Jsm;

        Dn.Si = Si;
        Dn.Ss = Ss;
        Dn.So = So;
        Dn.Sno = Sno;
        Dn.Snh = Snh;
        Dn.Snd = Snd;
        Dn.Salk = Salk;

        // upward connections
        Up.vS_dn = vS;
        Up.X_dn = X;
        annotation (
          Documentation(info="This class models the layers between the influent layer (feed_layer) and the lowest layer (bottom_layer)
of an ASM1 secondary clarifier based on Haertel.

Hydraulic and sedimentation flux (mass exchange in same direction) with above and underneath layer.

Sedimentation flux is calculated based on the sedimentation velocity
function and the omega correction function by Haertel.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end lower_layer;

      model feed_layer "Influent layer of Haertel`s SC model"
        import WWSC = WasteWater.BSM1.SecClar.Haertel.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;

        WWU.MassConcentration Xf "sludge concentration in clarifier feed";
        SI.Length z "vertical coordinate of current layer";

        parameter SI.Length hsc;
        parameter Integer i "number of layers above feed layer";
        Real omega;

        WWSC.LowerLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        WWSC.UpperLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
        BSM1.Interfaces.WWFlowAsm1in In
          annotation (Placement(transformation(extent={{-110,-6},{-90,14}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, ISV);
        omega = WWSC.omega(z, Xf, hsc, zm, ISV, i);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = if vS < Dn.vS_dn then omega*(vS*X) else omega*min(vS*X, Dn.vS_dn*Dn.
          X_dn);

        // ODE of solid component
        der(X) = (In.Q/Asc*Xf - (-Up.Qe)/Asc*X - (-(Dn.Qr + Dn.Qw))/Asc*X + Up.
          SedFlux - Jsm)/zm;

        // ODE of soluble components
        der(Si) = (In.Q*In.Si - (-Up.Qe)*Si - (-(Dn.Qr + Dn.Qw))*Si)/(Asc*zm);
        der(Ss) = (In.Q*In.Ss - (-Up.Qe)*Ss - (-(Dn.Qr + Dn.Qw))*Ss)/(Asc*zm);
        der(So) = (In.Q*In.So - (-Up.Qe)*So - (-(Dn.Qr + Dn.Qw))*So)/(Asc*zm);
        der(Sno) = (In.Q*In.Sno - (-Up.Qe)*Sno - (-(Dn.Qr + Dn.Qw))*Sno)/(Asc*zm);
        der(Snh) = (In.Q*In.Snh - (-Up.Qe)*Snh - (-(Dn.Qr + Dn.Qw))*Snh)/(Asc*zm);
        der(Snd) = (In.Q*In.Snd - (-Up.Qe)*Snd - (-(Dn.Qr + Dn.Qw))*Snd)/(Asc*zm);
        der(Salk) = (In.Q*In.Salk - (-Up.Qe)*Salk - (-(Dn.Qr + Dn.Qw))*Salk)/(Asc*
          zm);

        // volume flow rates
        In.Q + Up.Qe + Dn.Qr + Dn.Qw = 0;

        Dn.SedFlux = -Jsm;
        Dn.X = X;

        Dn.Si = Si;
        Dn.Ss = Ss;
        Dn.So = So;
        Dn.Sno = Sno;
        Dn.Snh = Snh;
        Dn.Snd = Snd;
        Dn.Salk = Salk;

        Up.X_dn = X;

        Up.Si = Si;
        Up.Ss = Ss;
        Up.So = So;
        Up.Sno = Sno;
        Up.Snh = Snh;
        Up.Snd = Snd;
        Up.Salk = Salk;

        annotation (
          Documentation(info="This class models the influent layer of an ASM1 secondary clarifier based on Haertel.

It receives the wastewater stream from the biological part (feed).
Hydraulic and sedimentation flux (mass exchange in opposite directions) with above layer
and hydraulic and sedimentation flux (mass exchange in same direction) with underneath layer.

Sedimentation flux is calculated based on the sedimentation velocity
function and the omega correction function by Haertel.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end feed_layer;

      model upper_layer "Layer above influent of Haertels`s SC"

        import WWSC = WasteWater.BSM1.SecClar.Haertel.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWSC.UpperLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        WWSC.UpperLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, ISV);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = vS*X;

        // ODE of solid component
        der(X) = (Dn.Qe/Asc*(Dn.X_dn - X) + Up.SedFlux - Jsm)/zm;

        // ODEs of soluble components
        der(Si) = Dn.Qe*(Dn.Si - Si)/(Asc*zm);
        der(Ss) = Dn.Qe*(Dn.Ss - Ss)/(Asc*zm);
        der(So) = Dn.Qe*(Dn.So - So)/(Asc*zm);
        der(Sno) = Dn.Qe*(Dn.Sno - Sno)/(Asc*zm);
        der(Snh) = Dn.Qe*(Dn.Snh - Snh)/(Asc*zm);
        der(Snd) = Dn.Qe*(Dn.Snd - Snd)/(Asc*zm);
        der(Salk) = Dn.Qe*(Dn.Salk - Salk)/(Asc*zm);

        // downward connection
        Dn.SedFlux = -Jsm;

        // upward connections
        Up.Qe + Dn.Qe = 0;

        Up.X_dn = X;

        Up.Si = Si;
        Up.Ss = Ss;
        Up.So = So;
        Up.Sno = Sno;
        Up.Snh = Snh;
        Up.Snd = Snd;
        Up.Salk = Salk;

        annotation (
          Documentation(info="This class models the layers between the influent layer (feed_layer) and the effluent layer (top_layer)
of an ASM1 secondary clarifier based on Haertel.

Hydraulic and sedimentation flux (mass exchange in opposite directions) with above and underneath layer.

Sedimentation flux is calculated based on the sedimentation velocity
function by Haertel."),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end upper_layer;

      model top_layer "Effluent layer of Haertel`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Haertel.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        extends WWSC.ratios;
        WWSC.UpperLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        BSM1.Interfaces.WWFlowAsm1out Out
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, ISV);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = vS*X;

        // ODE of solid component
        der(X) = (Dn.Qe/Asc*(Dn.X_dn - X) - Jsm)/zm;

        // ODEs of soluble components
        der(Si) = Dn.Qe*(Dn.Si - Si)/(Asc*zm);
        der(Ss) = Dn.Qe*(Dn.Ss - Ss)/(Asc*zm);
        der(So) = Dn.Qe*(Dn.So - So)/(Asc*zm);
        der(Sno) = Dn.Qe*(Dn.Sno - Sno)/(Asc*zm);
        der(Snh) = Dn.Qe*(Dn.Snh - Snh)/(Asc*zm);
        der(Snd) = Dn.Qe*(Dn.Snd - Snd)/(Asc*zm);
        der(Salk) = Dn.Qe*(Dn.Salk - Salk)/(Asc*zm);

        // downward connection
        Dn.SedFlux = -Jsm;

        // effluent volume flow rate
        Out.Q + Dn.Qe = 0;

        // effluent, solid and soluble components (ASM1)
        Out.Si = Si;
        Out.Ss = Ss;
        Out.Xi = rXi*X;
        Out.Xs = rXs*X;
        Out.Xbh = rXbh*X;
        Out.Xba = rXba*X;
        Out.Xp = rXp*X;
        Out.So = So;
        Out.Sno = Sno;
        Out.Snh = Snh;
        Out.Snd = Snd;
        Out.Xnd = rXnd*X;
        Out.Salk = Salk;

        annotation (
          Documentation(info="This class models the top layer of an ASM1 secondary clarifier based on Haertel.

No sedimentation flux (mass exchange) with above but hydraulic and sedimentation flux
(opposite directions) underneath.
From here effluent goes to the receiving water.

Sedimentation flux is calculated based on the sedimentation velocity
function by Haertel.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,58},{-8,40},{10,40},{10,32},{22,50},{10,66},{10,58},{-8,
                    58}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,58},{-8,40},{10,40},{10,32},{22,50},{10,66},{10,58},{-8,
                    58}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end top_layer;
      annotation (
        Documentation(info="This package contains classes (layer models) to built ASM1 secondary clarifier models, an Interfaces sub-library
and provides an ASM1 10-layer secondary clarifier model all bases on Haertel`s [1]
sedimentation velocity and omega correction functions.

A secondary clarifier layer model needs at least a top_layer, a feed_layer and a bottom_layer
and may have several upper_layer in between above the feed_layer and several lower_layer in
between below the feed_layer.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

References:

[1] L. Haertel: Modellansaetze zur dynamischen Simulation des Belebtschlammverfahrens.
    TH Darmstadt, Dissertation, 1990.

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2002 - 2003, Gerald Reichl
"));
    end Haertel;

    package Krebs "Secondary settling tank modelling by Krebs (ASM1)"
      extends Modelica.Icons.Library;

      package Interfaces
        "Partial models for Secondary Clarifier Model by Krebs"

        extends Modelica.Icons.Library;

        partial model SCVar "partial models providing variables"
          package WWU = WasteWater.WasteWaterUnits;
          WWU.MassConcentration Xf "total sludge concentration";
          WWU.MassConcentration XB "sludge concentration in sludge layer";
          WWU.MassConcentration XR "sludge concentration of return";

          WWU.MassConcentration Si1(fixed=true)
            "Soluble inert organic matter in first stirrer tank of the excess layer";
          WWU.MassConcentration Ss1(fixed=true)
            "Readily biodegradable substrate in first stirrer tank of the excess layer";
          WWU.MassConcentration So1(fixed=true)
            "Dissolved oxygen in first stirrer tank of the excess layer";
          WWU.MassConcentration Sno1(fixed=true)
            "Nitrate and nitrite nitrogen in first stirrer tank of the excess layer";
          WWU.MassConcentration Snh1(fixed=true)
            "Ammonium nitrogen in first stirrer tank of the excess layer";
          WWU.MassConcentration Snd1(fixed=true)
            "Soluble biodegradable organic nitrogen in first stirrer tank of the excess layer";
          WWU.Alkalinity Salk1(fixed=true)
            "Alkalinity in first stirrer tank of the excess layer";

          WWU.MassConcentration Si2(fixed=true)
            "Soluble inert organic matter in second stirrer tank of the excess layer";
          WWU.MassConcentration Ss2(fixed=true)
            "Readily biodegradable substrate in second stirrer tank of the excess layer";
          WWU.MassConcentration So2(fixed=true)
            "Dissolved oxygen in second stirrer tank of the excess layer";
          WWU.MassConcentration Sno2(fixed=true)
            "Nitrate and nitrite nitrogen in second stirrer tank of the excess layer";
          WWU.MassConcentration Snh2(fixed=true)
            "Ammonium nitrogen in second stirrer tank of the excess layer";
          WWU.MassConcentration Snd2(fixed=true)
            "Soluble biodegradable organic nitrogen in second stirrer tank of the excess layer";
          WWU.Alkalinity Salk2(fixed=true)
            "Alkalinity in second stirrer tank of the excess layer";
          annotation (
            Documentation(info="partial models providing ASM1 variables"));
        end SCVar;

        partial model ratios "partial model for ratios of solid components"
          // ratios of solid components
          Real rXi;
          Real rXs;
          Real rXbh;
          Real rXba;
          Real rXp;
          Real rXnd;
          annotation (
            Documentation(info="partial model for ASM1 ratios of solid components"));
        end ratios;
        annotation (
          Documentation(info="This package contains partial models for ASM1 secondary clarifier models.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2003, Gerald Reichl
"));
      end Interfaces;

      model SecClarModKrebs "ASM1 Secondary Settling Tank Model based on Krebs"

        extends WasteWater.Icons.SecClarKrebs;
        import WWSC = WasteWater.BSM1.SecClar.Krebs.Interfaces;
        extends WWSC.SCVar;
        extends WWSC.ratios;

        import SI = Modelica.SIunits;
        package WI = WasteWater.BSM1.Interfaces;
        package WWU = WasteWater.WasteWaterUnits;
        parameter SI.Length hsc=4.0 "height of secondary clarifier";
        parameter SI.Area Asc=1500.0 "area of secondary clarifier";
        parameter WWU.SludgeVolumeIndex ISV=130 "Sludge Volume Index";
        Real te "thickening time in sludge layer in [d]";
        SI.Length hs "height of sludge layer";
        SI.Length he "height of excess layer";
        WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},
                  {-90,24}})));
        WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,
                  47},{112,67}})));
        WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
                  -106},{-20,-86}})));
        WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,
                  -106},{40,-86}})));
      equation

        // total sludge concentration in clarifier feed
        Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

        // ratios of solid components
        rXs = Feed.Xs/Xf;
        rXbh = Feed.Xbh/Xf;
        rXba = Feed.Xba/Xf;
        rXp = Feed.Xp/Xf;
        rXi = Feed.Xi/Xf;
        rXnd = Feed.Xnd/Xf;

          //following expression is only for steady state initial equation of XB and is necessary

          //to calculate hs, if there would be problems with "initial()" in your modelica version
        //leave out this term and initial XB (or hs) manually e.g. via script-file
        if initial() then
          XB = Feed.Q/(0.7*(-(Return.Q + Waste.Q)))*Xf;
        end if;

        //thickening time in sludge layer in [d]
        te = 5/7*Asc*hs/(-(Return.Q + Waste.Q));

        //sludge concentration in sludge layer (unit of time in [h]) in [g/m3]
        XB = (1000/ISV*((te*24)^(1/3)))*1000;

        //sludge concentration of return
        XR = 0.7*XB;

        //ODE of height of sludge layer
        der(hs) = (Feed.Q*Xf - (-(Return.Q + Waste.Q))*XR)/(Asc/2*XB);

        //height of excess layer
        he = hsc - hs;

        // ODE of soluble components in first stirrer tank of the excess layer
        der(Si1) = (Feed.Q*Feed.Si - (-Effluent.Q)*Si1 - (-(Waste.Q + Return.Q))*
          Si1)/(Asc*he/2);
        der(Ss1) = (Feed.Q*Feed.Ss - (-Effluent.Q)*Ss1 - (-(Waste.Q + Return.Q))*
          Ss1)/(Asc*he/2);
        der(So1) = (Feed.Q*Feed.So - (-Effluent.Q)*So1 - (-(Waste.Q + Return.Q))*
          So1)/(Asc*he/2);
        der(Sno1) = (Feed.Q*Feed.Sno - (-Effluent.Q)*Sno1 - (-(Waste.Q + Return.Q))
          *Sno1)/(Asc*he/2);
        der(Snh1) = (Feed.Q*Feed.Snh - (-Effluent.Q)*Snh1 - (-(Waste.Q + Return.Q))
          *Snh1)/(Asc*he/2);
        der(Snd1) = (Feed.Q*Feed.Snd - (-Effluent.Q)*Snd1 - (-(Waste.Q + Return.Q))
          *Snd1)/(Asc*he/2);
        der(Salk1) = (Feed.Q*Feed.Salk - (-Effluent.Q)*Salk1 - (-(Waste.Q + Return.
          Q))*Salk1)/(Asc*he/2);

        // ODE of soluble components in second stirrer tank of the excess layer
        der(Si2) = ((-Effluent.Q)*Si1 - (-Effluent.Q)*Si2)/(Asc*he/2);
        der(Ss2) = ((-Effluent.Q)*Ss1 - (-Effluent.Q)*Ss2)/(Asc*he/2);
        der(So2) = ((-Effluent.Q)*So1 - (-Effluent.Q)*So2)/(Asc*he/2);
        der(Sno2) = ((-Effluent.Q)*Sno1 - (-Effluent.Q)*Sno2)/(Asc*he/2);
        der(Snh2) = ((-Effluent.Q)*Snh1 - (-Effluent.Q)*Snh2)/(Asc*he/2);
        der(Snd2) = ((-Effluent.Q)*Snd1 - (-Effluent.Q)*Snd2)/(Asc*he/2);
        der(Salk2) = ((-Effluent.Q)*Salk1 - (-Effluent.Q)*Salk2)/(Asc*he/2);

        // volume flow rates
        Feed.Q + Effluent.Q + Return.Q + Waste.Q = 0;

        // effluent, solid and soluble components (ASM1)
        Effluent.Si = Si2;
        Effluent.Ss = Ss2;
        Effluent.So = So2;
        Effluent.Sno = Sno2;
        Effluent.Snh = Snh2;
        Effluent.Snd = Snd2;
        Effluent.Salk = Salk2;
        Effluent.Xi = 0.0*XR;
        Effluent.Xs = 0.0*XR;
        Effluent.Xbh = 0.0*XR;
        Effluent.Xba = 0.0*XR;
        Effluent.Xp = 0.0*XR;
        Effluent.Xnd = 0.0*XR;

        // return sludge flow, solid and soluble components (ASM1)
        Return.Si = Si1;
        Return.Ss = Ss1;
        Return.So = So1;
        Return.Sno = Sno1;
        Return.Snh = Snh1;
        Return.Snd = Snd1;
        Return.Salk = Salk1;
        Return.Xi = rXi*XR;
        Return.Xs = rXs*XR;
        Return.Xbh = rXbh*XR;
        Return.Xba = rXba*XR;
        Return.Xp = rXp*XR;
        Return.Xnd = rXnd*XR;

        // waste sludge flow, solid and soluble components (ASM1)
        Waste.Si = Si1;
        Waste.Ss = Ss1;
        Waste.So = So1;
        Waste.Sno = Sno1;
        Waste.Snh = Snh1;
        Waste.Snd = Snd1;
        Waste.Salk = Salk1;
        Waste.Xi = rXi*XR;
        Waste.Xs = rXs*XR;
        Waste.Xbh = rXbh*XR;
        Waste.Xba = rXba*XR;
        Waste.Xp = rXp*XR;
        Waste.Xnd = rXnd*XR;
        annotation (
          Documentation(info="This component models an ASM1 secondary clarifier based on Krebs conceptional model.
It consists of two compartments: a \"sludge-bed\" and a clear water zone above.

Parameters:
  hsc -  height of clarifier [m]
  Asc -  surface area of secondary clarifier [m2]
  ISV -  Sludge Volume Index [ml/g]
"),       Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(
                extent={{-90,80},{92,14}},
                lineColor={0,0,255},
                lineThickness=0.5),
              Rectangle(
                extent={{-90,14},{92,-86}},
                lineColor={0,0,255},
                lineThickness=0.5),
              Polygon(
                points={{-8,-20},{-8,-38},{-16,-38},{0,-48},{16,-38},{8,-38},{8,-20},
                    {-8,-20}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,34},{-8,54},{-16,54},{0,64},{16,54},{8,54},{8,34},{-8,
                    34}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Text(extent={{-90,78},{-34,66}}, textString=
                                                   "top_layer"),
              Text(extent={{-90,20},{-30,-16}}, textString=
                                                    "bottom_layer"),
              Line(
                points={{-90,48},{92,48}},
                color={0,0,255},
                pattern=LinePattern.Dash)}));
      end SecClarModKrebs;
      annotation (
        Documentation(info="This package contains an ASM1 secondary clarifier model and an Interfaces sub-library
based on Krebs conceptional model [1].
The settler model consists of two compartments, a \"sludge-bed\" and a clear water zone above.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

References:

[1] P. Krebs and M. Armbruster and W. Rodi: Numerische Nachklaerbeckenmodelle. Korrespondenz Abwasser. 47 (7)
    2000. pp 985-999.

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2003, Gerald Reichl
"));
    end Krebs;

    package Otterpohl "Secondary settling tank modelling by Otterpohl"
      extends Modelica.Icons.Library;

      package Interfaces
        "Connectors and partial ASM1 models for Secondary Clarifier Model by Otterpohl"

        extends Modelica.Icons.Library;

        connector UpperLayerPin "Connector above influent layer"

          package WWU = WasteWater.WasteWaterUnits;
          // effluent flow
          flow WWU.VolumeFlowRate Qe;
          // sedimentation flux (from micro and macro flocs)
          flow WWU.SedimentationFlux SedFlux_F;
          // caused by macro flocs
          flow WWU.SedimentationFlux SedFlux_S;
          // caused by micro flocs

            // sludge concentration of macro and micro flocs in (m-1)-th layer (dn=down)
          WWU.MassConcentration X_dn_F;
          WWU.MassConcentration X_dn_S;

          // soluble components
          WWU.MassConcentration Si;
          WWU.MassConcentration Ss;
          WWU.MassConcentration So;
          WWU.MassConcentration Sno;
          WWU.MassConcentration Snh;
          WWU.MassConcentration Snd;
          WWU.Alkalinity Salk;
          annotation (
            Documentation(info=
                  "Connector for ASM1 information and mass exchange between layers above the influent layer (feed_layer)."));
        end UpperLayerPin;

        connector LowerLayerPin "Connector below influent layer"

          package WWU = WasteWater.WasteWaterUnits;

          // return and waste sludge flow Qr, Qw
          flow WWU.VolumeFlowRate Qr;
          flow WWU.VolumeFlowRate Qw;

          // sedimentation flux (from micro and macro flocs)
          flow WWU.SedimentationFlux SedFlux_F;
          // caused by macro flocs
          flow WWU.SedimentationFlux SedFlux_S;
          // caused by micro flocs

          // total sludge concentration of micro and macro flocs in m-th layer
          WWU.MassConcentration X_F;
          WWU.MassConcentration X_S;

            // total sludge concentration of micro and macro flocs in (m-1)-th layer (dn=down)
          WWU.MassConcentration X_dn_F;
          WWU.MassConcentration X_dn_S;
          // sink velocity of macro flocs in (m-1)-th layer
          WWU.SedimentationVelocity vS_dn_F;
          // soluble components
          WWU.MassConcentration Si;
          WWU.MassConcentration Ss;
          WWU.MassConcentration So;
          WWU.MassConcentration Sno;
          WWU.MassConcentration Snh;
          WWU.MassConcentration Snd;
          WWU.Alkalinity Salk;
          annotation (
            Documentation(info=
                  "Connector for ASM1 information and mass exchange between layers below the influent layer (feed_layer)."));
        end LowerLayerPin;

        partial model SCParam "partial model providing clarifier parameters"
          import SI = Modelica.SIunits;
          package WWU = WasteWater.WasteWaterUnits;
          parameter SI.Length zm;
          parameter SI.Area Asc;
          parameter WWU.SludgeVolumeIndex ISV;
          parameter WWU.SedimentationVelocity vS_S=0.24;
          // 0.01[m/h]*24 -> [m/d]

          annotation (
            Documentation(info="partial model providing clarifier parameters"));
        end SCParam;

        partial model SCVar "partial models providing variables"

          package WWU = WasteWater.WasteWaterUnits;
          WWU.MassConcentration X "total sludge concentration in m-th layer";
          WWU.MassConcentration X_F "sludge concentration of macro flocs";
          WWU.MassConcentration X_S "sludge concentration of micro flocs";
          WWU.SedimentationVelocity vS_F "sink velocity of makro flocs";
          WWU.SedimentationFlux Jsm_F "sedimentation flux of macro flocs";
          WWU.SedimentationFlux Jsm_S "sedimentation flux of micro flocs";

          WWU.MassConcentration Si "Soluble inert organic matter";
          WWU.MassConcentration Ss "Readily biodegradable substrate";
          WWU.MassConcentration So "Dissolved oxygen";
          WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
          WWU.MassConcentration Snh "Ammonium nitrogen";
          WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
          WWU.Alkalinity Salk "Alkalinity";
          annotation (
            Documentation(info="partial models providing ASM1 variables"));
        end SCVar;

        partial model ratios "partial model for ratios of solid components"

          // ratios of solid components
          Real rXi;
          Real rXs;
          Real rXbh;
          Real rXba;
          Real rXp;
          Real rXnd;
          annotation (
            Documentation(info="partial model for ASM1 ratios of solid components"));
        end ratios;

        function vSfun "Sedimentation velocity function"

          // total sludge concentration in m-th layer in g/m3 or mg/l
          input Real X;
          //Sludge Volume Index
          input Real ISV;
          // sink velocity in m/d
          output Real vS;
        protected
          Real v0 "maximum settling velocity";
          Real nv "exponent as part of the Vesilind equation";
        algorithm
          v0 := (17.4*(exp(-0.0113*ISV)) + 3.931)*24;
          //[m/d]
          nv := (-0.9834*(exp(-0.00581*ISV)) + 1.043);
          //[l/g]
          vS := v0*exp(-nv*X/1000);
          annotation (
            Documentation(info="Sedimentation velocity function"));
        end vSfun;

        function omega "Omega correction function by Haertel"

          input Real z;
          //vertical coordinate, bottom: z=0
          input Real Xf;
          // total sludge concentration in clarifier feed
          input Real hsc;
          //height of secondary clarifier
          input Real zm;
          //height of m-th secondary clarifier layer
          input Real ISV;
          //Sludge Volume Index
          input Integer i;
          //number of layers above feed layer

          // correction function omega by Haertel based on [g/l]
          output Real omega;

        protected
          Real Xc "solids concentration at compression point";
          Real nv "exponent as part of the Vesilind equation";
          Real ht "height of transition point";
          Real hc "height of compressing point";
          Real B3;
          Real B4;

        algorithm
          Xc := 480/ISV;
          nv := 1.043 - 0.9834*exp(-0.00581*ISV);
          hc := (Xf/1000)*(hsc - zm*(i + 0.5))/Xc*(1.0 - 1.0/(Xc*nv));
          // unit change
          ht := min(2.0*hc, hsc - zm*(i + 0.5));

          B4 := 1.0 + 2.0*ISV/(100.0 + ISV);
          B3 := -((2*ISV + 100.0)/ISV)*hc^B4;

          omega := (1.0 - B3*ht^(-B4))/(1.0 - B3*z^(-B4));
          omega := min(1.0, omega);
          annotation (
            Documentation(info=
                  "This is Haertels omega correction function for the settling process."));
        end omega;
        annotation (
          Documentation(info="This package contains connectors and interfaces (partial models) for
the ASM1 secondary clarifier model based on Otterpohl [1] (two settling velocities for
distinction between micro and macro flocs and omega correction function).

References:

[1] R. Otterpohl and M. Freund: Dynamic models for clarifiers of activated sludge plants
    with dry and wet weather flows. Water Science and Technology. 26 (1992), pp 1391-1400.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2003, Gerald Reichl
"));
      end Interfaces;

      model SecClarModOtter
        "Secondary Clarifier Model based on Otterpohl (ASM1)"
        import WasteWater;

        extends WasteWater.Icons.SecClar;
        extends Interfaces.ratios;
        package SCP = Otterpohl;
        import SI = Modelica.SIunits;
        package WI = WasteWater.BSM1.Interfaces;
        package WWU = WasteWater.WasteWaterUnits;
        parameter SI.Length hsc=4.0 "height of secondary clarifier";
        parameter Integer n=10 "number of layers of SC model";
        parameter SI.Length zm=hsc/(1.0*n)
          "height of m-th secondary clarifier layer";
        parameter SI.Area Asc=1500.0 "area of secondary clarifier";
        parameter WWU.SludgeVolumeIndex ISV=130 "Sludge Volume Index";
        parameter Integer i=2
          "number of layers above current feed layer in this model";

        // total sludge concentration in clarifier feed
        WWU.MassConcentration Xf;

        // layers 1 to 10
        SCP.bottom_layer S1(
          zm=zm,
          Asc=Asc,
          ISV=ISV,
          rXs=rXs,
          rXbh=rXbh,
          rXba=rXba,
          rXp=rXp,
          rXi=rXi,
          rXnd=rXnd) annotation (Placement(transformation(extent={{-35,-93},{35,-78}})));
        SCP.lower_layer S2(
          hsc=hsc,
          zm=zm,
          z=(zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-74},{35,-59}})));
        SCP.lower_layer S3(
          hsc=hsc,
          zm=zm,
          z=(2*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-55},{35,-40}})));
        SCP.lower_layer S4(
          hsc=hsc,
          zm=zm,
          z=(3*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-36},{35,-21}})));
        SCP.lower_layer S5(
          hsc=hsc,
          zm=zm,
          z=(4*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-17},{35,-2}})));
        SCP.lower_layer S6(
          hsc=hsc,
          zm=zm,
          z=(5*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,2},{35,17}})));
        SCP.lower_layer S7(
          hsc=hsc,
          zm=zm,
          z=(6*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,21},{35,36}})));
        SCP.feed_layer S8(
          hsc=hsc,
          zm=zm,
          z=(7*zm + zm/2),
          Asc=Asc,
          ISV=ISV,
          i=i,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,40},{35,55}})));
        SCP.upper_layer S9(
          zm=zm,
          Asc=Asc,
          ISV=ISV) annotation (Placement(transformation(extent={{-35,59},{35,74}})));
        SCP.top_layer S10(
          zm=zm,
          Asc=Asc,
          ISV=ISV,
          rXs=rXs,
          rXbh=rXbh,
          rXba=rXba,
          rXp=rXp,
          rXi=rXi,
          rXnd=rXnd) annotation (Placement(transformation(extent={{-35,78},{35,93}})));
        WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},
                  {-90,24}})));
        WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,
                  47},{112,67}})));
        WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
                  -106},{-20,-86}})));
        WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,
                  -106},{40,-86}})));
      equation

        connect(S1.Up, S2.Dn) annotation (Line(points={{-2.22045e-15,-78},{
                -2.22045e-15,-74}}));
        connect(S2.Up, S3.Dn) annotation (Line(points={{-2.22045e-15,-59},{
                -2.22045e-15,-55}}));
        connect(S3.Up, S4.Dn) annotation (Line(points={{-2.22045e-15,-40},{
                -2.22045e-15,-36}}));
        connect(S5.Up, S6.Dn) annotation (Line(points={{-2.22045e-15,-2},{
                -2.22045e-15,2}}));
        connect(S6.Up, S7.Dn) annotation (Line(points={{-2.22045e-15,17},{
                -2.22045e-15,21}}));
        connect(S7.Up, S8.Dn) annotation (Line(points={{-2.22045e-15,36},{
                -2.22045e-15,40}}));
        connect(S9.Up, S10.Dn) annotation (Line(points={{-2.22045e-15,74},{
                -2.22045e-15,78}}));
        connect(S4.Up, S5.Dn) annotation (Line(points={{-2.22045e-15,-21},{
                -2.22045e-15,-17}}));
        connect(S8.Up, S9.Dn) annotation (Line(points={{-2.22045e-15,55},{
                -2.22045e-15,59}}));
        connect(Feed, S8.In) annotation (Line(points={{-98,14},{-98,47.8},{-33,47.8}}));
        connect(S1.PQw, Waste) annotation (Line(points={{17.5,-93},{17.5,-100},{30,
                -100}}));
        connect(S10.Out, Effluent) annotation (Line(points={{35,85.5},{67.5,85.5},{
                67.5,57},{100,57}}));
        connect(S1.PQr, Return) annotation (Line(points={{-21,-93},{-21,-100},{-30,
                -100}}));

        // total sludge concentration in clarifier feed
        Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

        // ratios of solid components
        rXs = Feed.Xs/Xf;
        rXbh = Feed.Xbh/Xf;
        rXba = Feed.Xba/Xf;
        rXp = Feed.Xp/Xf;
        rXi = Feed.Xi/Xf;
        rXnd = Feed.Xnd/Xf;

        annotation (
          Documentation(info="This component models an ASM1 10 - layer secondary clarifier model with 4 layers above the feed_layer (including top_layer)
and 5 layers below the feed_layer (including bottom_layer) based on Otterpohl`s theory.

Parameters:
  hsc -  height of clarifier [m]
  n   -  number of layers
  Asc -  surface area of sec. clar. [m2]
  ISV -  Sludge Volume Index [ml/g]
  i   -  number of layers above feed layer
"));
      end SecClarModOtter;

      model bottom_layer "Bottom layer of Otterpohls`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Otterpohl.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        extends WWSC.ratios;
        BSM1.Interfaces.WWFlowAsm1out PQr
          annotation (Placement(transformation(extent={{-70,-110},{-50,-90}})));
        BSM1.Interfaces.WWFlowAsm1out PQw
          annotation (Placement(transformation(extent={{40,-110},{60,-90}})));
        WWSC.LowerLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
      equation

        // sink velocity
        vS_F = WWSC.vSfun(X_F, ISV);

        // sedimentation flux in bottom layer
        Jsm_F = 0.0;
        Jsm_S = 0.0;

        // ODE of solid component
        der(X_F) = ((Up.Qr + Up.Qw)/Asc*(Up.X_F - X_F) + Up.SedFlux_F)/zm;
        der(X_S) = ((Up.Qr + Up.Qw)/Asc*(Up.X_S - X_S) + Up.SedFlux_S)/zm;

        X = X_F + X_S;

        // ODEs of soluble components
        der(Si) = (Up.Qr + Up.Qw)*(Up.Si - Si)/(Asc*zm);
        der(Ss) = (Up.Qr + Up.Qw)*(Up.Ss - Ss)/(Asc*zm);
        der(So) = (Up.Qr + Up.Qw)*(Up.So - So)/(Asc*zm);
        der(Sno) = (Up.Qr + Up.Qw)*(Up.Sno - Sno)/(Asc*zm);
        der(Snh) = (Up.Qr + Up.Qw)*(Up.Snh - Snh)/(Asc*zm);
        der(Snd) = (Up.Qr + Up.Qw)*(Up.Snd - Snd)/(Asc*zm);
        der(Salk) = (Up.Qr + Up.Qw)*(Up.Salk - Salk)/(Asc*zm);

        // upward connection
        Up.vS_dn_F = vS_F;
        Up.X_dn_F = X_F;
        Up.X_dn_S = X_S;

        // return and waste sludge volume flow rates
        PQr.Q + Up.Qr = 0;
        PQw.Q + Up.Qw = 0;

        // return sludge flow, solid and soluble components (ASM1)
        PQr.Si = Si;
        PQr.Ss = Ss;
        PQr.Xi = rXi*X;
        PQr.Xs = rXs*X;
        PQr.Xbh = rXbh*X;
        PQr.Xba = rXba*X;
        PQr.Xp = rXp*X;
        PQr.So = So;
        PQr.Sno = Sno;
        PQr.Snh = Snh;
        PQr.Snd = Snd;
        PQr.Xnd = rXnd*X;
        PQr.Salk = Salk;

        // waste sludge flow, solid and soluble components (ASM1)
        PQw.Si = Si;
        PQw.Ss = Ss;
        PQw.Xi = rXi*X;
        PQw.Xs = rXs*X;
        PQw.Xbh = rXbh*X;
        PQw.Xba = rXba*X;
        PQw.Xp = rXp*X;
        PQw.So = So;
        PQw.Sno = Sno;
        PQw.Snh = Snh;
        PQw.Snd = Snd;
        PQw.Xnd = rXnd*X;
        PQw.Salk = Salk;

        annotation (
          Documentation(info="This class models the lowest layer of an ASM1 secondary clarifier based on Otterpohl.

No sedimentation flux (mass exchange) with underneath but hydraulic and sedimentation flux (same direction)
with above layer.
From here return and waste sludge is removed.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end bottom_layer;

      model lower_layer "Layer below influent of Otterpohl`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Otterpohl.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWU.MassConcentration Xf "sludge concentration in clarifier feed";
        SI.Length z "vertical coordinate of current layer";
        parameter SI.Length hsc;
        parameter Integer i "number of layers above feed layer";
        Real omega;
        WWSC.LowerLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
        WWSC.LowerLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
      equation

        // sink velocity
        vS_F = WWSC.vSfun(X_F, ISV);
        omega = WWSC.omega(z, Xf, hsc, zm, ISV, i);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm_F = if vS_F < Dn.vS_dn_F then omega*(vS_F*X_F) else omega*min(vS_F*X_F,
            Dn.vS_dn_F*Dn.X_dn_F);
        Jsm_S = omega*min(vS_S*X_S, vS_S*Dn.X_dn_S);

        // ODE of solid component
        der(X_F) = ((Up.Qr + Up.Qw)/Asc*(Up.X_F - X_F) + Up.SedFlux_F - Jsm_F)/zm;
        der(X_S) = ((Up.Qr + Up.Qw)/Asc*(Up.X_S - X_S) + Up.SedFlux_S - Jsm_S)/zm;

        X = X_F + X_S;

        // ODEs of soluble components
        der(Si) = (Up.Qr + Up.Qw)*(Up.Si - Si)/(Asc*zm);
        der(Ss) = (Up.Qr + Up.Qw)*(Up.Ss - Ss)/(Asc*zm);
        der(So) = (Up.Qr + Up.Qw)*(Up.So - So)/(Asc*zm);
        der(Sno) = (Up.Qr + Up.Qw)*(Up.Sno - Sno)/(Asc*zm);
        der(Snh) = (Up.Qr + Up.Qw)*(Up.Snh - Snh)/(Asc*zm);
        der(Snd) = (Up.Qr + Up.Qw)*(Up.Snd - Snd)/(Asc*zm);
        der(Salk) = (Up.Qr + Up.Qw)*(Up.Salk - Salk)/(Asc*zm);

        // downward connections
        Dn.Qr + Up.Qr = 0;
        Dn.Qw + Up.Qw = 0;

        Dn.X_F = X_F;
        Dn.X_S = X_S;
        Dn.SedFlux_F = -Jsm_F;
        Dn.SedFlux_S = -Jsm_S;

        Dn.Si = Si;
        Dn.Ss = Ss;
        Dn.So = So;
        Dn.Sno = Sno;
        Dn.Snh = Snh;
        Dn.Snd = Snd;
        Dn.Salk = Salk;

        // upward connections
        Up.vS_dn_F = vS_F;
        Up.X_dn_F = X_F;
        Up.X_dn_S = X_S;
        annotation (
          Documentation(info="This class models the layers between the influent layer (feed_layer) and the lowest layer (bottom_layer)
of an ASM1 secondary clarifier based on Otterpohl.

Hydraulic and sedimentation flux (mass exchange in same direction) with above and underneath layer.

Sedimentation flux is calculated based on two sedimentation velocities
(for macro and micro flocs) and the omega correction function by Haertel.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end lower_layer;

      model feed_layer "Influent layer of Otterpohl`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Otterpohl.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWU.MassConcentration Xf "sludge concentration in clarifier feed";
        SI.Length z "vertical coordinate of current layer";
        parameter SI.Length hsc;
        parameter Integer i "number of layers above feed layer";
        Real omega;
        Real fl;

        WWSC.LowerLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        WWSC.UpperLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
        BSM1.Interfaces.WWFlowAsm1in In
          annotation (Placement(transformation(extent={{-110,-6},{-90,14}})));
      equation

        // sink velocity
        vS_F = WWSC.vSfun(X_F, ISV);
        omega = WWSC.omega(z, Xf, hsc, zm, ISV, i);
        fl = (9.4/ISV)*exp(-1.07*Xf/1000);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm_F = if vS_F < Dn.vS_dn_F then omega*(vS_F*X_F) else omega*min(vS_F*X_F,
            Dn.vS_dn_F*Dn.X_dn_F);
        Jsm_S = omega*min(vS_S*X_S, vS_S*Dn.X_dn_S);

        // ODE of solid component
        der(X_F) = (In.Q/Asc*Xf*(1 - fl) - (-Up.Qe)/Asc*X_F - (-(Dn.Qr + Dn.Qw))/
          Asc*X_F + Up.SedFlux_F - Jsm_F)/zm;
        der(X_S) = (In.Q/Asc*Xf*fl - (-Up.Qe)/Asc*X_S - (-(Dn.Qr + Dn.Qw))/Asc*X_S
           + Up.SedFlux_S - Jsm_S)/zm;

        X = X_F + X_S;

        // ODE of soluble components
        der(Si) = (In.Q*In.Si - (-Up.Qe)*Si - (-(Dn.Qr + Dn.Qw))*Si)/(Asc*zm);
        der(Ss) = (In.Q*In.Ss - (-Up.Qe)*Ss - (-(Dn.Qr + Dn.Qw))*Ss)/(Asc*zm);
        der(So) = (In.Q*In.So - (-Up.Qe)*So - (-(Dn.Qr + Dn.Qw))*So)/(Asc*zm);
        der(Sno) = (In.Q*In.Sno - (-Up.Qe)*Sno - (-(Dn.Qr + Dn.Qw))*Sno)/(Asc*zm);
        der(Snh) = (In.Q*In.Snh - (-Up.Qe)*Snh - (-(Dn.Qr + Dn.Qw))*Snh)/(Asc*zm);
        der(Snd) = (In.Q*In.Snd - (-Up.Qe)*Snd - (-(Dn.Qr + Dn.Qw))*Snd)/(Asc*zm);
        der(Salk) = (In.Q*In.Salk - (-Up.Qe)*Salk - (-(Dn.Qr + Dn.Qw))*Salk)/(Asc*
          zm);

        // volume flow rates
        In.Q + Up.Qe + Dn.Qr + Dn.Qw = 0;

        Dn.SedFlux_F = -Jsm_F;
        Dn.SedFlux_S = -Jsm_S;
        Dn.X_F = X_F;
        Dn.X_S = X_S;

        Dn.Si = Si;
        Dn.Ss = Ss;
        Dn.So = So;
        Dn.Sno = Sno;
        Dn.Snh = Snh;
        Dn.Snd = Snd;
        Dn.Salk = Salk;

        Up.X_dn_F = X_F;
        Up.X_dn_S = X_S;

        Up.Si = Si;
        Up.Ss = Ss;
        Up.So = So;
        Up.Sno = Sno;
        Up.Snh = Snh;
        Up.Snd = Snd;
        Up.Salk = Salk;

        annotation (
          Documentation(info="This class models the influent layer of an ASM1 secondary clarifier based on Otterpohl.

It receives the wastewater stream from the biological part (feed).
Hydraulic and sedimentation flux (mass exchange in opposite directions) with above layer
and hydraulic and sedimentation flux (mass exchange in same direction) with underneath layer.

Sedimentation flux is calculated based on two sedimentation velocities
(for macro and micro flocs) and the omega correction function by Haertel.
"),       Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end feed_layer;

      model upper_layer "Layer above influent of Otterpohl`s SC"

        import WWSC = WasteWater.BSM1.SecClar.Otterpohl.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWSC.UpperLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        WWSC.UpperLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
      equation

        // sink velocity
        vS_F = WWSC.vSfun(X_F, ISV);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm_F = vS_F*X_F;
        Jsm_S = vS_S*X_S;

        // ODE of solid component
        der(X_F) = (Dn.Qe/Asc*(Dn.X_dn_F - X_F) + Up.SedFlux_F - Jsm_F)/zm;
        der(X_S) = (Dn.Qe/Asc*(Dn.X_dn_S - X_S) + Up.SedFlux_S - Jsm_S)/zm;

        X = X_F + X_S;

        // ODEs of soluble components
        der(Si) = Dn.Qe*(Dn.Si - Si)/(Asc*zm);
        der(Ss) = Dn.Qe*(Dn.Ss - Ss)/(Asc*zm);
        der(So) = Dn.Qe*(Dn.So - So)/(Asc*zm);
        der(Sno) = Dn.Qe*(Dn.Sno - Sno)/(Asc*zm);
        der(Snh) = Dn.Qe*(Dn.Snh - Snh)/(Asc*zm);
        der(Snd) = Dn.Qe*(Dn.Snd - Snd)/(Asc*zm);
        der(Salk) = Dn.Qe*(Dn.Salk - Salk)/(Asc*zm);

        // downward connection
        Dn.SedFlux_F = -Jsm_F;
        Dn.SedFlux_S = -Jsm_S;

        // upward connections
        Up.Qe + Dn.Qe = 0;

        Up.X_dn_F = X_F;
        Up.X_dn_S = X_S;

        Up.Si = Si;
        Up.Ss = Ss;
        Up.So = So;
        Up.Sno = Sno;
        Up.Snh = Snh;
        Up.Snd = Snd;
        Up.Salk = Salk;

        annotation (
          Documentation(info="This class models the layers between the influent layer (feed_layer) and the effluent layer (top_layer)
of an ASM1 secondary clarifier based on Otterpohl.

Hydraulic and sedimentation flux (mass exchange in opposite directions) with above and underneath layer.

Sedimentation flux is calculated based on two sedimentation velocities
(for macro and micro flocs)."),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end upper_layer;

      model top_layer "Effluent layer of Otterpohl`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Otterpohl.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        extends WWSC.ratios;
        WWSC.UpperLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        BSM1.Interfaces.WWFlowAsm1out Out
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      equation

        // sink velocity
        vS_F = WWSC.vSfun(X_F, ISV);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm_F = vS_F*X_F;
        Jsm_S = vS_S*X_S;

        // ODE of solid component
        der(X_F) = (Dn.Qe/Asc*(Dn.X_dn_F - X_F) - Jsm_F)/zm;
        der(X_S) = (Dn.Qe/Asc*(Dn.X_dn_S - X_S) - Jsm_S)/zm;

        X = X_F + X_S;

        // ODEs of soluble components
        der(Si) = Dn.Qe*(Dn.Si - Si)/(Asc*zm);
        der(Ss) = Dn.Qe*(Dn.Ss - Ss)/(Asc*zm);
        der(So) = Dn.Qe*(Dn.So - So)/(Asc*zm);
        der(Sno) = Dn.Qe*(Dn.Sno - Sno)/(Asc*zm);
        der(Snh) = Dn.Qe*(Dn.Snh - Snh)/(Asc*zm);
        der(Snd) = Dn.Qe*(Dn.Snd - Snd)/(Asc*zm);
        der(Salk) = Dn.Qe*(Dn.Salk - Salk)/(Asc*zm);

        // downward connection
        Dn.SedFlux_F = -Jsm_F;
        Dn.SedFlux_S = -Jsm_S;

        // effluent volume flow rate
        Out.Q + Dn.Qe = 0;

        // effluent, solid and soluble components (ASM1)
        Out.Si = Si;
        Out.Ss = Ss;
        Out.Xi = rXi*X;
        Out.Xs = rXs*X;
        Out.Xbh = rXbh*X;
        Out.Xba = rXba*X;
        Out.Xp = rXp*X;
        Out.So = So;
        Out.Sno = Sno;
        Out.Snh = Snh;
        Out.Snd = Snd;
        Out.Xnd = rXnd*X;
        Out.Salk = Salk;
        annotation (
          Documentation(info="This class models the top layer of an ASM1 secondary clarifier based on Otterpohl.

No sedimentation flux (mass exchange) with above but hydraulic and sedimentation flux
(opposite directions) underneath.
From here effluent goes to the receiving water.

Sedimentation flux is calculated based on two sedimentation velocities
(for micro and macro flocs).
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,58},{-8,40},{10,40},{10,32},{22,50},{10,66},{10,58},{-8,
                    58}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,58},{-8,40},{10,40},{10,32},{22,50},{10,66},{10,58},{-8,
                    58}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end top_layer;
      annotation (
        Documentation(info="This package contains classes (layer models) to built ASM1 secondary clarifier models, an Interfaces sub-library
and provides an ASM1 10-layer secondary clarifier model all bases on Otterpohls`s [1]
sedimentation velocities for macro and micro flocs and the omega correction function.

A secondary clarifier layer model needs at least a top_layer, a feed_layer and a bottom_layer
and may have several upper_layer in between above the feed_layer and several lower_layer in
between below the feed_layer.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

References:

[1] R. Otterpohl and M. Freund: Dynamic models for clarifiers of activated sludge plants
    with dry and wet weather flows. Water Science and Technology. 26 (1992), pp 1391-1400.

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2003, Gerald Reichl
"));
    end Otterpohl;

    package Simple "Simple ASM1 Secondary clarifier model"
      extends Modelica.Icons.Library;

      model SimpleSecClarMod "Simple ASM1 Secondary Clarifier Model"

        extends WasteWater.Icons.SecClarSimple;
        extends WasteWater.BSM1.SecClar.Takacs.Interfaces.ratios;
        import SI = Modelica.SIunits;
        package WI = WasteWater.BSM1.Interfaces;
        package WWU = WasteWater.WasteWaterUnits;

        parameter SI.Length hsc=4.0 "height of secondary clarifier";
        parameter SI.Area Asc=1500.0 "area of secondary clarifier";

        WWU.MassConcentration Xf "total sludge concentration in clarifier feed";
        WWU.MassConcentration X "sludge concentration in clarifier";
        WWU.MassConcentration Si "Soluble inert organic matter";
        WWU.MassConcentration Ss "Readily biodegradable substrate";
        WWU.MassConcentration So "Dissolved oxygen";
        WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
        WWU.MassConcentration Snh "Ammonium nitrogen";
        WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
        WWU.Alkalinity Salk "Alkalinity";

        WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},
                  {-90,24}})));
        WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,
                  47},{112,67}})));
        WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
                  -106},{-20,-86}})));
        WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,
                  -106},{40,-86}})));
      equation

        // total sludge concentration in clarifier feed
        Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

        // ratios of solid components
        rXs = Feed.Xs/Xf;
        rXbh = Feed.Xbh/Xf;
        rXba = Feed.Xba/Xf;
        rXp = Feed.Xp/Xf;
        rXi = Feed.Xi/Xf;
        rXnd = Feed.Xnd/Xf;

        // ODE of sludge concentration
        der(X) = (Feed.Q*Xf - (-(Waste.Q + Return.Q))*X)/(Asc*hsc);

        // ODE of soluble components
        der(Si) = (Feed.Q*Feed.Si - (-Effluent.Q)*Si - (-(Waste.Q + Return.Q))*Si)/
          (Asc*hsc);
        der(Ss) = (Feed.Q*Feed.Ss - (-Effluent.Q)*Ss - (-(Waste.Q + Return.Q))*Ss)/
          (Asc*hsc);
        der(So) = (Feed.Q*Feed.So - (-Effluent.Q)*So - (-(Waste.Q + Return.Q))*So)/
          (Asc*hsc);
        der(Sno) = (Feed.Q*Feed.Sno - (-Effluent.Q)*Sno - (-(Waste.Q + Return.Q))*
          Sno)/(Asc*hsc);
        der(Snh) = (Feed.Q*Feed.Snh - (-Effluent.Q)*Snh - (-(Waste.Q + Return.Q))*
          Snh)/(Asc*hsc);
        der(Snd) = (Feed.Q*Feed.Snd - (-Effluent.Q)*Snd - (-(Waste.Q + Return.Q))*
          Snd)/(Asc*hsc);
        der(Salk) = (Feed.Q*Feed.Salk - (-Effluent.Q)*Salk - (-(Waste.Q + Return.Q))
           *Salk)/(Asc*hsc);

        // volume flow rates
        Feed.Q + Effluent.Q + Return.Q + Waste.Q = 0;

        // effluent, solid and soluble components (ASM1)
        Effluent.Si = Si;
        Effluent.Ss = Ss;
        Effluent.Xi = 0.0*X;
        Effluent.Xs = 0.0*X;
        Effluent.Xbh = 0.0*X;
        Effluent.Xba = 0.0*X;
        Effluent.Xp = 0.0*X;
        Effluent.So = So;
        Effluent.Sno = Sno;
        Effluent.Snh = Snh;
        Effluent.Snd = Snd;
        Effluent.Xnd = 0.0*X;
        Effluent.Salk = Salk;

        // return sludge flow, solid and soluble components (ASM1)
        Return.Si = Si;
        Return.Ss = Ss;
        Return.Xi = rXi*X;
        Return.Xs = rXs*X;
        Return.Xbh = rXbh*X;
        Return.Xba = rXba*X;
        Return.Xp = rXp*X;
        Return.So = So;
        Return.Sno = Sno;
        Return.Snh = Snh;
        Return.Snd = Snd;
        Return.Xnd = rXnd*X;
        Return.Salk = Salk;

        // waste sludge flow, solid and soluble components (ASM1)
        Waste.Si = Si;
        Waste.Ss = Ss;
        Waste.Xi = rXi*X;
        Waste.Xs = rXs*X;
        Waste.Xbh = rXbh*X;
        Waste.Xba = rXba*X;
        Waste.Xp = rXp*X;
        Waste.So = So;
        Waste.Sno = Sno;
        Waste.Snh = Snh;
        Waste.Snd = Snd;
        Waste.Xnd = rXnd*X;
        Waste.Salk = Salk;

        annotation (
          Documentation(info="This component models very simple the secondary clarification process by
just using a single fully mixed tank which removes all particulate substances from the effluent
and returns the sludge. No sedimentation and compression, etc. is considered (for ASM1).

Parameters:
  hsc -    height of clarifier [m]
  Asc -    surface area of sec. clar. [m2]
"),       Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Polygon(
                points={{-20,-70},{20,-70},{4,-84},{-4,-84},{-20,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-4,-84},{4,-92}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-80,-48},{-36,-64},{38,-64},{80,-48},{-80,-48}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-80,62},{80,-40}},
                lineColor={223,191,159},
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Text(extent={{-80,98},{80,66}}, textString=
                                                  "%name"),
              Polygon(
                points={{-36,-64},{38,-64},{20,-70},{-20,-70},{-36,-64}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Line(
                points={{4,-92},{4,-84},{20,-70},{80,-48}},
                thickness=0.5),
              Rectangle(
                extent={{-80,-40},{80,-48}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{80,62},{92,54}},
                lineColor={0,127,255},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{80,54},{92,54}},
                thickness=0.5),
              Line(
                points={{-4,-92},{-4,-84},{-20,-70},{-80,-48},{-80,10}},
                thickness=0.5),
              Line(
                points={{-80,62},{-80,16}},
                thickness=0.5),
              Line(
                points={{-80,10},{-90,10}},
                thickness=0.5),
              Line(
                points={{-80,16},{-90,16}},
                thickness=0.5),
              Rectangle(
                extent={{-20,-92},{20,-98}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Line(
                points={{-20,-92},{-4,-92}},
                thickness=0.5),
              Line(
                points={{-20,-98},{20,-98}},
                thickness=0.5),
              Line(
                points={{20,-92},{4,-92}},
                thickness=0.5),
              Line(
                points={{80,-48},{80,54}},
                thickness=0.5),
              Text(extent={{-100,-60},{-40,-80}}, textString=
                                                      "return"),
              Text(extent={{40,-60},{100,-80}}, textString=
                                                    "waste"),
              Polygon(
                points={{16,44},{33,44},{31,52},{48,42},{31,31},{33,39},{16,39},{16,
                    44}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-46,32},{-29,32},{-31,40},{-14,30},{-31,19},{-29,27},{-46,
                    27},{-46,32}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{18,-26},{22,-26},{22,-42},{28,-40},{20,-54},{12,-40},{18,
                    -42},{18,-26}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={191,95,0},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-32,-10},{-28,-10},{-28,-26},{-22,-24},{-30,-38},{-38,-24},
                    {-32,-26},{-32,-10}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={191,95,0},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-90,16},{-80,10}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end SimpleSecClarMod;
      annotation (
        Documentation(info="This package just provides a very simple ASM1 secondary clarifier model
with no sludge storage, no sludge sedimentation and no use of layers.
The model consists of one tank removing all particulate substances.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2002, Gerald Reichl
"));
    end Simple;

    package Takacs "Secondary settling tank modelling by Takacs"
      extends Modelica.Icons.Library;

      package Interfaces
        "Connectors and partial models for Secondary Clarifier Model by Takacs"

        extends Modelica.Icons.Library;

        connector UpperLayerPin "Connector above influent layer"

          package WWU = WasteWater.WasteWaterUnits;

          // effluent flow
          flow WWU.VolumeFlowRate Qe;
          // sedimentation flux
          flow WWU.SedimentationFlux SedFlux;

            // total sludge concentration and sink velocity in (m-1)-th layer (dn=down)
          WWU.MassConcentration X_dn;
          WWU.SedimentationVelocity vS_dn;

          // soluble components
          WWU.MassConcentration Si;
          WWU.MassConcentration Ss;
          WWU.MassConcentration So;
          WWU.MassConcentration Sno;
          WWU.MassConcentration Snh;
          WWU.MassConcentration Snd;
          WWU.Alkalinity Salk;
          annotation (
            Documentation(info=
                  "Connector for ASM1 information and mass exchange between layers above the influent layer (feed_layer)."));
        end UpperLayerPin;

        connector LowerLayerPin "Connector below influent layer"

          package WWU = WasteWater.WasteWaterUnits;

          // return and waste sludge flow Qr, Qw
          flow WWU.VolumeFlowRate Qr;
          flow WWU.VolumeFlowRate Qw;

          // sedimentation flux
          flow WWU.SedimentationFlux SedFlux;

          // total sludge concentration in m-th layer
          WWU.MassConcentration X;

            // total sludge concentration and sink velocity in(m-1)-th layer (dn=down)
          WWU.MassConcentration X_dn;
          WWU.SedimentationVelocity vS_dn;

          // soluble components
          WWU.MassConcentration Si;
          WWU.MassConcentration Ss;
          WWU.MassConcentration So;
          WWU.MassConcentration Sno;
          WWU.MassConcentration Snh;
          WWU.MassConcentration Snd;
          WWU.Alkalinity Salk;
          annotation (
            Documentation(info=
                  "Connector for ASM1 information and mass exchange between layers below the influent layer (feed_layer)."));
        end LowerLayerPin;

        partial model SCParam "partial model providing clarifier parameters"
          import SI = Modelica.SIunits;
          parameter SI.Length zm;
          parameter SI.Area Asc;

          annotation (
            Documentation(info="partial model providing clarifier parameters"));
        end SCParam;

        partial model SCVar "partial models providing variables"

          package WWU = WasteWater.WasteWaterUnits;

          WWU.MassConcentration X "total sludge concentration in m-th layer";
          WWU.MassConcentration Xf
            "total sludge concentration in clarifier feed";
          WWU.SedimentationVelocity vS "sink velocity in m-th layer";
          WWU.SedimentationFlux Jsm "sedimentation flux m-th layer";

          WWU.MassConcentration Si "Soluble inert organic matter";
          WWU.MassConcentration Ss "Readily biodegradable substrate";
          WWU.MassConcentration So "Dissolved oxygen";
          WWU.MassConcentration Sno "Nitrate and nitrite nitrogen";
          WWU.MassConcentration Snh "Ammonium nitrogen";
          WWU.MassConcentration Snd "Soluble biodegradable organic nitrogen";
          WWU.Alkalinity Salk "Alkalinity";
          annotation (
            Documentation(info="partial models providing ASM1 variables"));
        end SCVar;

        partial model ratios "partial model for ratios of solid components"

          // ratios of solid components
          Real rXi;
          Real rXs;
          Real rXbh;
          Real rXba;
          Real rXp;
          Real rXnd;
          annotation (
            Documentation(info="partial model for ASM1 ratios of solid components"));
        end ratios;

        function vSfun "Sedimentation velocity function"

          // total sludge concentration in m-th layer in g/m3 or mg/l
          input Real X;
          // total sludge concentration in clarifier feed in g/m3 or mg/l
          input Real Xf;
          // sink velocity in m/d
          output Real vS;

          parameter Real v0slash=250.0 "max. settling velocity in m/d";
          parameter Real v0=474.0 "max. Vesilind settl. veloc. in m/d";
          parameter Real rh=0.000576 "hindered zone settl. param. in m3/(g SS)";
          parameter Real rp=0.00286
            "flocculant zone settl. param. in m3/(g SS)";
          parameter Real fns=0.00228 "non-settleable fraction in -";
        algorithm

          // computation of sink velocity
          vS := max(0.0, min(v0slash, v0*(exp(-rh*(X - fns*Xf)) - exp(-rp*(X - fns*
            Xf)))));

          annotation (
            Documentation(info=
                  "Takacs double-exponential sedimentation velocity function."));
        end vSfun;
        annotation (
          Documentation(info="This package contains connectors and interfaces (partial models) for
the ASM1 secondary clarifier model based on Takacs [1] (double-exponential settling velocity).

References:

[1] I. Takacs and G.G. Patry and D. Nolasco: A dynamic model of the clarification-thickening
    process. Water Research. 25 (1991) 10, pp 1263-1271.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2001, Gerald Reichl
"));
      end Interfaces;

      model SecClarModTakacs "Secondary Clarifier ASM1 Model based on Takacs"

        extends WasteWater.Icons.SecClar;
        extends Interfaces.ratios;
        package SCP = Takacs;
        import SI = Modelica.SIunits;
        package WI = WasteWater.BSM1.Interfaces;
        package WWU = WasteWater.WasteWaterUnits;

        parameter SI.Length hsc=4.0 "height of secondary clarifier";
        parameter Integer n=10 "number of layers of SC model";
        parameter SI.Length zm=hsc/(1.0*n)
          "height of m-th secondary clarifier layer";
        parameter SI.Area Asc=1500.0 "area of secondary clarifier";
        parameter WWU.MassConcentration Xt=3000.0 "threshold for X";

        // total sludge concentration in clarifier feed
        WWU.MassConcentration Xf;

        WI.WWFlowAsm1in Feed annotation (Placement(transformation(extent={{-110,4},
                  {-90,24}})));
        WI.WWFlowAsm1out Effluent annotation (Placement(transformation(extent={{92,
                  47},{112,67}})));
        WI.WWFlowAsm1out Return annotation (Placement(transformation(extent={{-40,
                  -106},{-20,-86}})));
        WI.WWFlowAsm1out Waste annotation (Placement(transformation(extent={{20,
                  -106},{40,-86}})));

        // layers 1 to 10
        SCP.bottom_layer S1(
          zm=zm,
          Asc=Asc,
          Xf=Xf,
          rXs=rXs,
          rXbh=rXbh,
          rXba=rXba,
          rXp=rXp,
          rXi=rXi,
          rXnd=rXnd) annotation (Placement(transformation(extent={{-35,-93},{35,-78}})));
        SCP.lower_layer S2(
          zm=zm,
          Asc=Asc,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-74},{35,-59}})));
        SCP.lower_layer S3(
          zm=zm,
          Asc=Asc,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-55},{35,-40}})));
        SCP.lower_layer S4(
          zm=zm,
          Asc=Asc,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-36},{35,-21}})));
        SCP.lower_layer S5(
          zm=zm,
          Asc=Asc,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,-17},{35,-2}})));
        SCP.feed_layer S6(
          zm=zm,
          Asc=Asc,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,2},{35,17}})));
        SCP.upper_layer S7(
          zm=zm,
          Asc=Asc,
          Xf=Xf,
          Xt=Xt) annotation (Placement(transformation(extent={{-35,21},{35,36}})));
        SCP.upper_layer S8(
          zm=zm,
          Asc=Asc,
          Xt=Xt,
          Xf=Xf) annotation (Placement(transformation(extent={{-35,40},{35,55}})));
        SCP.upper_layer S9(
          zm=zm,
          Asc=Asc,
          Xf=Xf,
          Xt=Xt) annotation (Placement(transformation(extent={{-35,59},{35,74}})));
        SCP.top_layer S10(
          zm=zm,
          Asc=Asc,
          Xf=Xf,
          Xt=Xt,
          rXs=rXs,
          rXbh=rXbh,
          rXba=rXba,
          rXp=rXp,
          rXi=rXi,
          rXnd=rXnd) annotation (Placement(transformation(extent={{-35,78},{35,93}})));
      equation

        connect(S1.Up, S2.Dn) annotation (Line(points={{-2.22045e-15,-78},{
                -2.22045e-15,-74}}));
        connect(S2.Up, S3.Dn) annotation (Line(points={{-2.22045e-15,-59},{
                -2.22045e-15,-55}}));
        connect(S3.Up, S4.Dn) annotation (Line(points={{-2.22045e-15,-40},{
                -2.22045e-15,-36}}));
        connect(S5.Up, S6.Dn) annotation (Line(points={{-2.22045e-15,-2},{
                -2.22045e-15,2}}));
        connect(S6.Up, S7.Dn) annotation (Line(points={{-2.22045e-15,17},{
                -2.22045e-15,21}}));
        connect(S7.Up, S8.Dn) annotation (Line(points={{-2.22045e-15,36},{
                -2.22045e-15,40}}));
        connect(S9.Up, S10.Dn) annotation (Line(points={{-2.22045e-15,74},{
                -2.22045e-15,78}}));
        connect(S4.Up, S5.Dn) annotation (Line(points={{-2.22045e-15,-21},{
                -2.22045e-15,-17}}));
        connect(S8.Up, S9.Dn) annotation (Line(points={{-2.22045e-15,55},{
                -2.22045e-15,59}}));
        connect(Feed, S6.In) annotation (Line(points={{-100,10},{-67.5,10},{-67.5,
                9.8},{-35,9.8}}));
        connect(S1.PQw, Waste) annotation (Line(points={{17.5,-93},{17.5,-100},{30,
                -100}}));
        connect(S10.Out, Effluent) annotation (Line(points={{35,85.5},{67.5,85.5},{
                67.5,57},{100,57}}));
        connect(S1.PQr, Return) annotation (Line(points={{-21,-93},{-21,-100},{-30,
                -100}}));

        // total sludge concentration in clarifier feed
        Xf = 0.75*(Feed.Xs + Feed.Xbh + Feed.Xba + Feed.Xp + Feed.Xi);

        // ratios of solid components
        rXs = Feed.Xs/Xf;
        rXbh = Feed.Xbh/Xf;
        rXba = Feed.Xba/Xf;
        rXp = Feed.Xp/Xf;
        rXi = Feed.Xi/Xf;
        rXnd = Feed.Xnd/Xf;
        annotation (
          Documentation(info="This component models an ASM1 10 - layer secondary clarifier model with 4 layers above the feed_layer (including top_layer)
and 5 layers below the feed_layer (including bottom_layer) based on Takac`s theory.

Parameters:
  hsc -  height of clarifier [m]
  n   -  number of layers
  Asc -  surface area of sec. clar. [m2]
  Xt  -  threshold value for Xtss [mg/l]

"));
      end SecClarModTakacs;

      model bottom_layer "Bottom layer of Takac`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Takacs.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        extends WWSC.ratios;
        BSM1.Interfaces.WWFlowAsm1out PQr
          annotation (Placement(transformation(extent={{-70,-110},{-50,-90}})));
        BSM1.Interfaces.WWFlowAsm1out PQw
          annotation (Placement(transformation(extent={{40,-110},{60,-90}})));
        WWSC.LowerLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, Xf);
        Jsm = 0.0;

        // ODE of solid component
        der(X) = ((Up.Qr + Up.Qw)/Asc*(Up.X - X) + Up.SedFlux)/zm;

        // ODEs of soluble components
        der(Si) = (Up.Qr + Up.Qw)*(Up.Si - Si)/(Asc*zm);
        der(Ss) = (Up.Qr + Up.Qw)*(Up.Ss - Ss)/(Asc*zm);
        der(So) = (Up.Qr + Up.Qw)*(Up.So - So)/(Asc*zm);
        der(Sno) = (Up.Qr + Up.Qw)*(Up.Sno - Sno)/(Asc*zm);
        der(Snh) = (Up.Qr + Up.Qw)*(Up.Snh - Snh)/(Asc*zm);
        der(Snd) = (Up.Qr + Up.Qw)*(Up.Snd - Snd)/(Asc*zm);
        der(Salk) = (Up.Qr + Up.Qw)*(Up.Salk - Salk)/(Asc*zm);

        // upward connection
        Up.vS_dn = vS;
        Up.X_dn = X;

        // return and waste sludge volume flow rates
        PQr.Q + Up.Qr = 0;
        PQw.Q + Up.Qw = 0;

        // return sludge flow, solid and soluble components (ASM1)
        PQr.Si = Si;
        PQr.Ss = Ss;
        PQr.Xi = rXi*X;
        PQr.Xs = rXs*X;
        PQr.Xbh = rXbh*X;
        PQr.Xba = rXba*X;
        PQr.Xp = rXp*X;
        PQr.So = So;
        PQr.Sno = Sno;
        PQr.Snh = Snh;
        PQr.Snd = Snd;
        PQr.Xnd = rXnd*X;
        PQr.Salk = Salk;

        // waste sludge flow, solid and soluble components (ASM1)
        PQw.Si = Si;
        PQw.Ss = Ss;
        PQw.Xi = rXi*X;
        PQw.Xs = rXs*X;
        PQw.Xbh = rXbh*X;
        PQw.Xba = rXba*X;
        PQw.Xp = rXp*X;
        PQw.So = So;
        PQw.Sno = Sno;
        PQw.Snh = Snh;
        PQw.Snd = Snd;
        PQw.Xnd = rXnd*X;
        PQw.Salk = Salk;

        annotation (
          Documentation(info="This class models the lowest layer of an ASM1 secondary clarifier based on Takacs.

No sedimentation flux (mass exchange) with underneath but hydraulic and sedimentation flux (same direction)
with above layer.
From here return and waste sludge is removed.
"),       Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end bottom_layer;

      model lower_layer "Layer below influent of Takac`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Takacs.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWSC.LowerLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
        WWSC.LowerLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, Xf);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = min(vS*X, Dn.vS_dn*Dn.X_dn);

        // ODE of solid component
        der(X) = ((Up.Qr + Up.Qw)/Asc*(Up.X - X) + Up.SedFlux - Jsm)/zm;

        // ODEs of soluble components
        der(Si) = (Up.Qr + Up.Qw)*(Up.Si - Si)/(Asc*zm);
        der(Ss) = (Up.Qr + Up.Qw)*(Up.Ss - Ss)/(Asc*zm);
        der(So) = (Up.Qr + Up.Qw)*(Up.So - So)/(Asc*zm);
        der(Sno) = (Up.Qr + Up.Qw)*(Up.Sno - Sno)/(Asc*zm);
        der(Snh) = (Up.Qr + Up.Qw)*(Up.Snh - Snh)/(Asc*zm);
        der(Snd) = (Up.Qr + Up.Qw)*(Up.Snd - Snd)/(Asc*zm);
        der(Salk) = (Up.Qr + Up.Qw)*(Up.Salk - Salk)/(Asc*zm);

        // downward connections
        Dn.Qr + Up.Qr = 0;
        Dn.Qw + Up.Qw = 0;

        Dn.X = X;
        Dn.SedFlux = -Jsm;

        Dn.Si = Si;
        Dn.Ss = Ss;
        Dn.So = So;
        Dn.Sno = Sno;
        Dn.Snh = Snh;
        Dn.Snd = Snd;
        Dn.Salk = Salk;

        // upward connections
        Up.vS_dn = vS;
        Up.X_dn = X;
        annotation (
          Documentation(info="This class models the layers between the influent layer (feed_layer) and the lowest layer (bottom_layer)
of an ASM1 secondary clarifier based on Takacs.

Hydraulic and sedimentation flux (mass exchange in same direction) with above and underneath layer.

Sedimentation flux is calculated based on the double-exponential sedimentation velocity
function by Takacs."),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,68},{-68,50},{-76,50},{-60,40},{-44,50},{-52,50},{-52,
                    68},{-68,68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end lower_layer;

      model feed_layer "Influent layer of Takac`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Takacs.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        WWSC.LowerLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        WWSC.UpperLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
        BSM1.Interfaces.WWFlowAsm1in In
          annotation (Placement(transformation(extent={{-110,-6},{-90,14}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, Xf);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = min(vS*X, Dn.vS_dn*Dn.X_dn);

        // ODE of solid component
        der(X) = (In.Q/Asc*Xf - (-Up.Qe)/Asc*X - (-(Dn.Qr + Dn.Qw))/Asc*X + Up.
          SedFlux - Jsm)/zm;

        // ODE of soluble components
        der(Si) = (In.Q*In.Si - (-Up.Qe)*Si - (-(Dn.Qr + Dn.Qw))*Si)/(Asc*zm);
        der(Ss) = (In.Q*In.Ss - (-Up.Qe)*Ss - (-(Dn.Qr + Dn.Qw))*Ss)/(Asc*zm);
        der(So) = (In.Q*In.So - (-Up.Qe)*So - (-(Dn.Qr + Dn.Qw))*So)/(Asc*zm);
        der(Sno) = (In.Q*In.Sno - (-Up.Qe)*Sno - (-(Dn.Qr + Dn.Qw))*Sno)/(Asc*zm);
        der(Snh) = (In.Q*In.Snh - (-Up.Qe)*Snh - (-(Dn.Qr + Dn.Qw))*Snh)/(Asc*zm);
        der(Snd) = (In.Q*In.Snd - (-Up.Qe)*Snd - (-(Dn.Qr + Dn.Qw))*Snd)/(Asc*zm);
        der(Salk) = (In.Q*In.Salk - (-Up.Qe)*Salk - (-(Dn.Qr + Dn.Qw))*Salk)/(Asc*
          zm);

        // volume flow rates
        In.Q + Up.Qe + Dn.Qr + Dn.Qw = 0;

        Dn.SedFlux = -Jsm;
        Dn.X = X;

        Dn.Si = Si;
        Dn.Ss = Ss;
        Dn.So = So;
        Dn.Sno = Sno;
        Dn.Snh = Snh;
        Dn.Snd = Snd;
        Dn.Salk = Salk;

        Up.vS_dn = vS;
        Up.X_dn = X;

        Up.Si = Si;
        Up.Ss = Ss;
        Up.So = So;
        Up.Sno = Sno;
        Up.Snh = Snh;
        Up.Snd = Snd;
        Up.Salk = Salk;
        annotation (
          Documentation(info="This class models the influent layer of an ASM1 secondary clarifier based on Takacs.

It receives the wastewater stream from the biological part (feed).
Hydraulic and sedimentation flux (mass exchange in opposite directions) with above layer
and hydraulic and sedimentation flux (mass exchange in same direction) with underneath layer.

Sedimentation flux is calculated based on the double-exponential sedimentation velocity
function by Takacs."),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-40},{-68,-58},{-76,-58},{-60,-68},{-44,-58},{-52,-58},
                    {-52,-40},{-68,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid)}));
      end feed_layer;

      model upper_layer "Layer above influent of Takac`s SC"
        // Xt = Xthreshold

        import WWSC = WasteWater.BSM1.SecClar.Takacs.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        parameter WWU.MassConcentration Xt;

        WWSC.UpperLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        WWSC.UpperLayerPin Up annotation (Placement(transformation(extent={{-10,90},
                  {10,110}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, Xf);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = if Dn.X_dn <= Xt then vS*X else min(vS*X, Dn.vS_dn*Dn.X_dn);

        // ODE of solid component
        der(X) = (Dn.Qe/Asc*(Dn.X_dn - X) + Up.SedFlux - Jsm)/zm;

        // ODEs of soluble components
        der(Si) = Dn.Qe*(Dn.Si - Si)/(Asc*zm);
        der(Ss) = Dn.Qe*(Dn.Ss - Ss)/(Asc*zm);
        der(So) = Dn.Qe*(Dn.So - So)/(Asc*zm);
        der(Sno) = Dn.Qe*(Dn.Sno - Sno)/(Asc*zm);
        der(Snh) = Dn.Qe*(Dn.Snh - Snh)/(Asc*zm);
        der(Snd) = Dn.Qe*(Dn.Snd - Snd)/(Asc*zm);
        der(Salk) = Dn.Qe*(Dn.Salk - Salk)/(Asc*zm);

        // downward connection
        Dn.SedFlux = -Jsm;

        // upward connections
        Up.Qe + Dn.Qe = 0;

        Up.vS_dn = vS;
        Up.X_dn = X;

        Up.Si = Si;
        Up.Ss = Ss;
        Up.So = So;
        Up.Sno = Sno;
        Up.Snh = Snh;
        Up.Snd = Snd;
        Up.Salk = Salk;
        annotation (
          Documentation(info="This class models the layers between the influent layer (feed_layer) and the effluent layer (top_layer)
an ASM1 secondary clarifier based on Takacs.

Hydraulic and sedimentation flux (mass exchange in opposite directions) with above and underneath layer.

Sedimentation flux is calculated based on the double-exponential sedimentation velocity
function by Takacs."),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,40},{-68,60},{-76,60},{-60,70},{-44,60},{-52,60},{-52,
                    40},{-68,40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{52,68},{52,50},{44,50},{60,40},{76,50},{68,50},{68,68},{52,
                    68}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end upper_layer;

      model top_layer "Effluent layer of Takac`s SC model"

        import WWSC = WasteWater.BSM1.SecClar.Takacs.Interfaces;
        extends WWSC.SCParam;
        extends WWSC.SCVar;
        extends WWSC.ratios;

        parameter WWU.MassConcentration Xt;
        // Xt = Xthreshold

        WWSC.UpperLayerPin Dn annotation (Placement(transformation(extent={{-10,
                  -110},{10,-90}})));
        BSM1.Interfaces.WWFlowAsm1out Out
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      equation

        // sink velocity
        vS = WWSC.vSfun(X, Xf);

        // sedimentation flux in m-th layer sinking to lower layer
        Jsm = if Dn.X_dn <= Xt then vS*X else min(vS*X, Dn.vS_dn*Dn.X_dn);

        // ODE of solid component
        der(X) = (Dn.Qe/Asc*(Dn.X_dn - X) - Jsm)/zm;

        // ODEs of soluble components
        der(Si) = Dn.Qe*(Dn.Si - Si)/(Asc*zm);
        der(Ss) = Dn.Qe*(Dn.Ss - Ss)/(Asc*zm);
        der(So) = Dn.Qe*(Dn.So - So)/(Asc*zm);
        der(Sno) = Dn.Qe*(Dn.Sno - Sno)/(Asc*zm);
        der(Snh) = Dn.Qe*(Dn.Snh - Snh)/(Asc*zm);
        der(Snd) = Dn.Qe*(Dn.Snd - Snd)/(Asc*zm);
        der(Salk) = Dn.Qe*(Dn.Salk - Salk)/(Asc*zm);

        // downward connection
        Dn.SedFlux = -Jsm;

        // effluent volume flow rate
        Out.Q + Dn.Qe = 0;

        // effluent, solid and soluble components (ASM1)
        Out.Si = Si;
        Out.Ss = Ss;
        Out.Xi = rXi*X;
        Out.Xs = rXs*X;
        Out.Xbh = rXbh*X;
        Out.Xba = rXba*X;
        Out.Xp = rXp*X;
        Out.So = So;
        Out.Sno = Sno;
        Out.Snh = Snh;
        Out.Snd = Snd;
        Out.Xnd = rXnd*X;
        Out.Salk = Salk;
        annotation (
          Documentation(info="This class models the top layer of an ASM1 secondary clarifier based on Takacs.

No sedimentation flux (mass exchange) with above but hydraulic and sedimentation flux
(opposite directions) underneath.
From here effluent goes to the receiving water.

Sedimentation flux is calculated based on the double-exponential sedimentation velocity
function by Takacs."),
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Text(extent={{-100,20},{100,-20}}, textString=
                                                     "%name"),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,58},{-8,40},{10,40},{10,32},{22,50},{10,66},{10,58},{-8,
                    58}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={2,2}), graphics={
              Rectangle(extent={{-100,100},{100,-100}}),
              Polygon(
                points={{52,-40},{52,-58},{44,-58},{60,-68},{76,-58},{68,-58},{68,
                    -40},{52,-40}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={223,191,159},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-8,58},{-8,40},{10,40},{10,32},{22,50},{10,66},{10,58},{-8,
                    58}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-68,-70},{-68,-50},{-76,-50},{-60,-40},{-44,-50},{-52,-50},
                    {-52,-70},{-68,-70}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end top_layer;
      annotation (
        Documentation(info="This package contains classes (layer models) to built ASM1 secondary clarifier models, an Interfaces sub-library
and provides an ASM1 10-layer secondary clarifier model all bases on Takacs [1]
double exponential sedimentation velocity function.

A secondary clarifier layer model needs at least a top_layer, a feed_layer and a bottom_layer
and may have several upper_layer in between above the feed_layer and several lower_layer in
between below the feed_layer.

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

References:

[1] I. Takacs and G.G. Patry and D. Nolasco: A dynamic model of the clarification-thickening
    process. Water Research. 25 (1991) 10, pp 1263-1271.

This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2002, Gerald Reichl
"));
    end Takacs;
  annotation (
    Documentation(info="This library provides a collection of ASM1 secondary clarifier models based on
several theories.

The library currently is structured in following sub-libraries.

 - Takacs: secondary clarifier model by Takacs et al [1]
 - Haertel: secondary clarifier model by Haertel [2]
 - Otterpohl: secondary clarifier model by Otterpohl [3]
 - Krebs: secondary clarifier model by Krebs [4]
 - Simple: very basic secondary clarifier model


Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de

References:

[1] I. Takacs and G.G. Patry and D. Nolasco: A dynamic model of the clarification-thickening
    process. Water Research. 25 (1991) 10, pp 1263-1271.

[2] L. Haertel: Modellansaetze zur dynamischen Simulation des Belebtschlammverfahrens.
    TH Darmstadt, Dissertation, 1990.

[3] R. Otterpohl and M. Freund: Dynamic models for clarifiers of activated sludge plants
    with dry and wet weather flows. Water Science and Technology. 26 (1992), pp 1391-1400.

[4] P. Krebs and M. Armbruster and W. Rodi: Numerische Nachklaerbeckenmodelle. Korrespondenz Abwasser. 47 (7)
    2000. pp 985-999.


This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2003, Gerald Reichl
"));
  end SecClar;

model nitri5 "ASM1 nitrification tank"
  // nitrification (aerated) tank, based on the ASM1 model

  extends WasteWater.Icons.nitri;
  extends Interfaces.ASM1base;

  // tank specific parameters
  parameter Modelica.SIunits.Volume V=1333 "Volume of nitrification tank";

  // aeration system dependent parameters
  parameter Real alpha=0.7 "Oxygen transfer factor";
  parameter Modelica.SIunits.Length de=4.5 "depth of aeration";
  parameter Real R_air=23.5 "specific oxygen feed factor [gO2/(m^3*m)]";
  WWU.MassConcentration So_sat "Dissolved oxygen saturation";
  //parameter Real Kla = 84;

  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-110,
            -10},{-90,10}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{90,
            -10},{110,10}})));
  Interfaces.WWFlowAsm1out MeasurePort annotation (Placement(transformation(
          extent={{50,40},{60,50}})));

  Modelica.Blocks.Interfaces.RealInput T annotation (Placement(transformation(
          extent={{-110,30},{-90,50}})));
  Modelica.Blocks.Interfaces.RealOutput Kla "Oxygen transfer factor" annotation (Placement(transformation(extent={{-10,-10},
              {10,10}},
          rotation=90,
          origin={-30,50})));
  Interfaces.AirFlow AirIn annotation (Placement(transformation(extent={{-5,
            -103},{5,-93}})));
equation

  // Temperature dependent oxygen saturation by Simba
  // So_sat =13.89 + (-0.3825 + (0.007311 - 0.00006588*T)*T)*T;

  So_sat = 8;
  Kla = 84;

  // extends the Oxygen differential equation by an aeration term
  // aeration [mgO2/l]; AirIn.Q_air needs to be in
  // Simulationtimeunit [m3*day^-1]
  // aeration = (alpha*(So_sat - So)/So_sat*AirIn.Q_air*R_air*de)/V;
  aeration = Kla * (So_sat - So);
  // So = 2;
  // volume dependent dilution term of each concentration

  inputSi = (In.Si - Si)*In.Q/V;
  inputSs = (In.Ss - Ss)*In.Q/V;
  inputXi = (In.Xi - Xi)*In.Q/V;
  inputXs = (In.Xs - Xs)*In.Q/V;
  inputXbh = (In.Xbh - Xbh)*In.Q/V;
  inputXba = (In.Xba - Xba)*In.Q/V;
  inputXp = (In.Xp - Xp)*In.Q/V;
  inputSo = (In.So - So)*In.Q/V;
  inputSno = (In.Sno - Sno)*In.Q/V;
  inputSnh = (In.Snh - Snh)*In.Q/V;
  inputSnd = (In.Snd - Snd)*In.Q/V;
  inputXnd = (In.Xnd - Xnd)*In.Q/V;
  inputSalk = (In.Salk - Salk)*In.Q/V;

  annotation (
    Documentation(info="This component models the ASM1 processes and reactions taking place in an aerated (nitrification) tank
of a wastewater treatment plant.

The InPort signal gives the tank temperature to the model.

Parameters:

        - all soichiometric and kinetic parameters of the activated sludge model No.1 (ASM1)
  V     - volume of the tank [m3]
  alpha - oxygen transfer factor
  de    - depth of the aeration system [m]
  R_air - specific oxygen feed factor [g O2/(m3*m)]
"), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}),
        graphics),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics));
end nitri5;

model sensor_aeration_energy
    "Ideal sensor to measure chemical oxygen demand (COD)"

  extends WasteWater.Icons.sensor_COD;

  Modelica.Blocks.Interfaces.RealInput Kla5 annotation (Placement(
        transformation(extent={{-20,-20},{20,20}},
          rotation=90,
          origin={0,-80}),                            iconTransformation(extent={{-16,-90},
              {4,-70}})));
  Modelica.Blocks.Interfaces.RealOutput AE(start=0) annotation (Placement(
        transformation(extent={{80,-10},{100,10}}), iconTransformation(extent={{
            80,-10},{100,10}})));
  Real T(start=1e-3);
  Real Kla3 = 240;
  Real Kla4 = 240;
  Real So_sat = 8;

equation
  der(T) = 1.0;
  der(AE*T) = So_sat/1.8/1000*1333*(Kla3 + Kla4 + Kla5);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics), Diagram(coordinateSystem(
            preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
          graphics));
end sensor_aeration_energy;

model sensor_effluent_quality "Influent quality index"
  extends WasteWater.Icons.sensor_COD;
  extends Interfaces.stoichiometry;

  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-80,-10},
            {-60,10}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput EQ( start=0) annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));

  Real T(start=1e-3);
  Real B_SS=2;
  Real B_COD=1;
  Real B_NKj=30;
  Real B_NO=10;
  Real B_BOD5=2;
  Real So_sat=8;
  Real S_NKj0;
  Real SS_0;
  Real BOD_50;
  Real COD_0;

equation
  In.Q + Out.Q = 0;

  In.Si = Out.Si;
  In.Ss = Out.Ss;
  In.Xi = Out.Xi;
  In.Xs = Out.Xs;
  In.Xbh = Out.Xbh;
  In.Xba = Out.Xba;
  In.Xp = Out.Xp;
  In.So = Out.So;
  In.Sno = Out.Sno;
  In.Snh = Out.Snh;
  In.Snd = Out.Snd;
  In.Xnd = Out.Xnd;
  In.Salk = Out.Salk;

  S_NKj0 = In.Snh + In.Snd + In.Xnd + i_xb * (In.Xbh + In.Xba) + i_xp * (In.Xp + In.Xi);
  SS_0 = 0.75 * (In.Xs + In.Xi + In.Xbh + In.Xba + In.Xp);
  BOD_50 = 0.65 * (In.Ss + In.Xs + (1-f_p) * (In.Xbh + In.Xba));
  COD_0 = In.Ss + In.Si + In.Xi + In.Xbh + In.Xba + In.Xp;
  der(T) = 1.0;
  der(EQ*T) = (1/1000) * (B_SS * SS_0 * B_COD * COD_0 + B_NKj*S_NKj0 + B_NO*In.Sno + B_BOD5*BOD_50)*In.Q;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end sensor_effluent_quality;

model sensor_pump_energy "Total pump energy"

  extends WasteWater.Icons.sensor_COD;

  Interfaces.WWFlowAsm1in In_a annotation (Placement(transformation(extent={{-70,50},
              {-50,70}})));
  Interfaces.WWFlowAsm1out Out_a annotation (Placement(transformation(extent={{-70,-70},
              {-50,-50}})));
  Interfaces.WWFlowAsm1in In_r annotation (Placement(transformation(extent={{-10,50},
              {10,70}})));
  Interfaces.WWFlowAsm1out Out_r annotation (Placement(transformation(extent={{-10,-70},
              {10,-50}})));
  Interfaces.WWFlowAsm1in In_w annotation (Placement(transformation(extent={{50,50},
              {70,70}})));
  Interfaces.WWFlowAsm1out Out_w annotation (Placement(transformation(extent={{50,-70},
              {70,-50}})));

  Modelica.Blocks.Interfaces.RealOutput PE( start = 0) annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));

  Real T(start=1e-3);
equation

  In_a.Q + Out_a.Q = 0;
  In_r.Q + Out_r.Q = 0;
  In_w.Q + Out_w.Q = 0;

  der(T) = 1.0;
  der(PE*T) = 0.004*In_a.Q + 0.008*In_r.Q + 0.05*In_w.Q;

  // eventually abs(In.Q) to be shure to have pos. signal
  In_a.Si = Out_a.Si;
  In_a.Ss = Out_a.Ss;
  In_a.Xi = Out_a.Xi;
  In_a.Xs = Out_a.Xs;
  In_a.Xbh = Out_a.Xbh;
  In_a.Xba = Out_a.Xba;
  In_a.Xp = Out_a.Xp;
  In_a.So = Out_a.So;
  In_a.Sno = Out_a.Sno;
  In_a.Snh = Out_a.Snh;
  In_a.Snd = Out_a.Snd;
  In_a.Xnd = Out_a.Xnd;
  In_a.Salk = Out_a.Salk;

  In_r.Si = Out_r.Si;
  In_r.Ss = Out_r.Ss;
  In_r.Xi = Out_r.Xi;
  In_r.Xs = Out_r.Xs;
  In_r.Xbh = Out_r.Xbh;
  In_r.Xba = Out_r.Xba;
  In_r.Xp = Out_r.Xp;
  In_r.So = Out_r.So;
  In_r.Sno = Out_r.Sno;
  In_r.Snh = Out_r.Snh;
  In_r.Snd = Out_r.Snd;
  In_r.Xnd = Out_r.Xnd;
  In_r.Salk = Out_r.Salk;

  In_w.Si = Out_w.Si;
  In_w.Ss = Out_w.Ss;
  In_w.Xi = Out_w.Xi;
  In_w.Xs = Out_w.Xs;
  In_w.Xbh = Out_w.Xbh;
  In_w.Xba = Out_w.Xba;
  In_w.Xp = Out_w.Xp;
  In_w.So = Out_w.So;
  In_w.Sno = Out_w.Sno;
  In_w.Snh = Out_w.Snh;
  In_w.Snd = Out_w.Snd;
  In_w.Xnd = Out_w.Xnd;
  In_w.Salk = Out_w.Salk;

    annotation (Placement(transformation(
        origin={0,-98},
        extent={{-10,-10},{10,10}},
        rotation=270)), Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics),
    Documentation(info="This component measures the flow of an ASM1 wastewater stream and provides
the result as output signal (to be further processed with blocks of
the Modelica.Blocks library).
"));
end sensor_pump_energy;

model sensor_influent_quality "Influent quality index"

  extends WasteWater.Icons.sensor_COD;
  extends Interfaces.stoichiometry;

  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-80,-10},
            {-60,10}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{-10,
              -100},{10,-80}})));
  Modelica.Blocks.Interfaces.RealOutput IQ( start=0) annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));

  Real T(start=1e-3);
  Real B_SS=2;
  Real B_COD=1;
  Real B_NKj=30;
  Real B_NO=10;
  Real B_BOD5=2;
  Real So_sat=8;
  Real S_NKj0;
  Real SS_0;
  Real BOD_50;
  Real COD_0;

equation
  In.Q + Out.Q = 0;

  In.Si = Out.Si;
  In.Ss = Out.Ss;
  In.Xi = Out.Xi;
  In.Xs = Out.Xs;
  In.Xbh = Out.Xbh;
  In.Xba = Out.Xba;
  In.Xp = Out.Xp;
  In.So = Out.So;
  In.Sno = Out.Sno;
  In.Snh = Out.Snh;
  In.Snd = Out.Snd;
  In.Xnd = Out.Xnd;
  In.Salk = Out.Salk;

  S_NKj0 = In.Snh + In.Snd + In.Xnd + i_xb * (In.Xbh + In.Xba) + i_xp * (In.Xp + In.Xi);
  SS_0 = 0.75 * (In.Xs + In.Xi + In.Xbh + In.Xba + In.Xp);
  BOD_50 = 0.65 * (In.Ss + In.Xs + (1-f_p) * (In.Xbh + In.Xba));
  COD_0 = In.Ss + In.Si + In.Xi + In.Xbh + In.Xba + In.Xp;
  der(T) = 1.0;
  der(IQ*T) = (1/1000) * (B_SS * SS_0 * B_COD * COD_0 + B_NKj*S_NKj0 + B_NO*In.Sno + B_BOD5*BOD_50)*In.Q;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end sensor_influent_quality;

model sensor_sludge_production "Sludge production"

  extends WasteWater.Icons.sensor_COD;
  extends Interfaces.stoichiometry;

  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-80,-10},
            {-60,10}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput SP( start=0) annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));

  Real T(start=1e-3);
equation

  In.Q + Out.Q = 0;

  In.Si = Out.Si;
  In.Ss = Out.Ss;
  In.Xi = Out.Xi;
  In.Xs = Out.Xs;
  In.Xbh = Out.Xbh;
  In.Xba = Out.Xba;
  In.Xp = Out.Xp;
  In.So = Out.So;
  In.Sno = Out.Sno;
  In.Snh = Out.Snh;
  In.Snd = Out.Snd;
  In.Xnd = Out.Xnd;
  In.Salk = Out.Salk;

  der(T) = 1.0;
  der(SP*T) = 0.75*(In.Xs + In.Xi + In.Xbh + In.Xba + In.Xp)*In.Q;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end sensor_sludge_production;

model sensor_mixing_energy "Influent quality index"

  extends WasteWater.Icons.sensor_COD;
  extends Interfaces.stoichiometry;

  Interfaces.WWFlowAsm1in In annotation (Placement(transformation(extent={{-80,-10},
            {-60,10}})));
  Interfaces.WWFlowAsm1out Out annotation (Placement(transformation(extent={{-10,
            -110},{10,-90}})));
  Modelica.Blocks.Interfaces.RealOutput ME( start=0) annotation (Placement(
        transformation(extent={{88,-10},{108,10}})));

  Real T(start=1e-3);
  Real V_3 = 1333;
  Real V_4 = 1333;
  Real V_5 = 1333;

equation
  der(T) = 1.0;
  if Kla < 20 then
    der(ME*T) = 24 * 0.005*(V_3 + V_4 + V_5);
  else
    der(ME*T) = 0;
  end if;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end sensor_mixing_energy;

model OCI "Sludge production"

  extends WasteWater.Icons.sensor_COD;
  extends Interfaces.stoichiometry;

  Modelica.Blocks.Interfaces.RealOutput OCI( start=0) annotation (Placement(
        transformation(extent={{70,-10},{90,10}})));
  Modelica.Blocks.Interfaces.RealInput PE annotation (Placement(
        transformation(extent={{-90,-30},{-70,-10}})));
  Modelica.Blocks.Interfaces.RealInput AE annotation (Placement(
        transformation(extent={{-90,50},{-70,70}})));
  Modelica.Blocks.Interfaces.RealInput SP annotation (Placement(
        transformation(extent={{-90,-70},{-70,-50}})));
  Modelica.Blocks.Interfaces.RealInput EC annotation (Placement(
        transformation(extent={{-90,10},{-70,30}})));
  //Real T(start=1e-3);
equation

  OCI = AE + PE + 5*SP + 3*EC;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end OCI;
annotation (
  Documentation(info="This library contains components to build models of biological municipal
wastewater treatment plants based on the Activated Sludge Model No.1 (ASM1) by the
International Association on Water Quality (IAWQ) [1,2].


The library currently is structured in following sub-libraries.

 - Interfaces (partial ASM1 models and connectors)
 - PreClar (primary clarifier models)
 - SecClar (several secondary settling tank models)
 - Examples (wastewater treatment plant models)

Main Author:
   Gerald Reichl
   Technische Universitaet Ilmenau
   Faculty of Informatics and Automation
   Department Dynamics and Simulation of ecological Systems
   P.O. Box 10 05 65
   98684 Ilmenau
   Germany
   email: gerald.reichl@tu-ilmenau.de


References:

[1] M. Henze and C.P.L. Grady Jr and W. Gujer and G.v.R. Marais and T. Matsuo:
    Activated Sludge Model No.1. IAWQ, 1987.
[2] M. Henze and W.Gujer and T. Mino and. M.v. Loosdrecht: Activated Sludge
    Models ASM1, ASM2, ASM2d, and ASM3. IWA Task Group on Mathematical Modelling
    for Design and Operation of Biological Wastewater Treatment, 2000.


This package is free software; it can be redistributed and/or modified under the terms of the Modelica license, see the license conditions and the accompanying
disclaimer in the documentation of package Modelica in file \"Modelica/package.mo\".

Copyright (C) 2000 - 2002, Gerald Reichl
"));
end BSM1;
