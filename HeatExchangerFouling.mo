package HeatExchangerFouling5 "Demo of a heat exchanger model"
  extends Modelica.Icons.ExamplesPackage;

  model HeatExchangerFoulingSimulation "simulation for the heat exchanger model"
  Modelica.SIunits.Temp_K T1in "Cold fluid inlet temperature";
  Modelica.SIunits.Temp_K T1out "Cold fluid outlet temperature";
  Modelica.SIunits.Temp_K T2in "Hot fluid inlet temperature";
  Modelica.SIunits.Temp_K T2out "Hot fluid outlet temperature";
  Modelica.SIunits.TemperatureDifference Delta1;
  Modelica.SIunits.TemperatureDifference Delta2;
  Modelica.SIunits.TemperatureDifference Delta;
  Modelica.SIunits.TemperatureDifference DeltaLmF;
  Modelica.SIunits.Temp_K Tf "film temperature";
  Modelica.SIunits.TemperatureDifference Tf2 "film temperature";
  Modelica.SIunits.QualityFactor F "Correction factor";
  Modelica.SIunits.Area At "Heat Transfer area";
  Modelica.SIunits.CoefficientOfHeatTransfer Uf "U subject to fouling";
  Modelica.SIunits.CoefficientOfHeatTransfer Ud;
  Modelica.SIunits.HeatFlowRate Qf "Heat flow subject to fouling";
  Modelica.SIunits.HeatFlowRate Qref;
  Modelica.SIunits.HeatFlowRate Qloss"Rate loss";
  Modelica.SIunits.HeatFlowRate Eloss"Energy loss";
  Modelica.SIunits.HeatFlowRate VF"$$";
  Modelica.SIunits.ReynoldsNumber Re_t "Reynolds number tube side";
  Modelica.SIunits.QualityFactor a;
  Modelica.SIunits.QualityFactor b;
  Modelica.SIunits.QualityFactor ga "gamma";
  Modelica.SIunits.ChemicalPotential E "Activation energy";
  Modelica.SIunits.MolarInternalEnergy R;
  Modelica.SIunits.ThermalInsulance Rf "thermal resistance due to fouling";
  Modelica.SIunits.ThermalInsulance Rf2 "thermal resistance due to fouling";
  Modelica.SIunits.ThermalInsulance Ra "thermal resistance due to fouling";
  Modelica.SIunits.QualityFactor e;
  Modelica.SIunits.QualityFactor tao;
  extends Modelica.Icons.Example;
  //Temperatures are not used, they are replaced by functions;
  //replaceable package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater;
  replaceable package Medium = Modelica.Media.Water.StandardWaterOnePhase;
  //package Medium = Modelica.Media.Incompressible.Examples.Essotherm650;
  
  HeatExchangerFouling5.BaseClasses.BasicHX HEX(
      
      redeclare package Medium_1 =
          Medium,
      redeclare package Medium_2 =
          Medium,
      redeclare model HeatTransfer_1 =
          Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.LocalPipeFlowHeatTransfer
          (alpha0=1000),
      redeclare model HeatTransfer_2 =
          Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer
          (alpha0=2000),
      T_start_1=304,
      T_start_2=300,
      Twall_start=300,
      area_h_1=4044.0693,
      area_h_2= 674.0124842,c_wall=500,
      crossArea_1= 4.5e-4,
      crossArea_2= 4.5e-4,
      dT=10,
      k_wall=100,
      length=20,
      m_flow_start_1=0.2,
      m_flow_start_2=0.2,
      modelStructure_1=Modelica.Fluid.Types.ModelStructure.av_b,
      modelStructure_2=Modelica.Fluid.Types.ModelStructure.a_vb,
      nNodes=20,
      perimeter_1=0.075,
      perimeter_2=0.075,
      rho_wall=900,
      s_wall=0.005,
      use_T_start= false) annotation (Placement(visible = true, transformation(extent = {{-30, -30}, {30, 30}}, rotation = 0)));
  
    Modelica.Fluid.Sources.Boundary_pT End2(nPorts= 1,
      p=1e5,
      T=280,
      redeclare package Medium = Medium) annotation (Placement(
          visible = true, transformation(origin = {-50, -60},extent = {{10, -10}, {-10, 10}}, rotation = 180)));
    Modelica.Fluid.Sources.Boundary_pT End1(nPorts= 1,
      p=1e5,
      T=300,
      redeclare package Medium = Medium) annotation (Placement(
          visible = true, transformation(extent = {{96, -10}, {76, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T Hot2(
      redeclare package Medium = Medium,
      T= 360,
      m_flow= 0.2,nPorts= 1,
      use_T_in= false,
      use_X_in=false,
      use_m_flow_in=true)
                  annotation (Placement(visible = true, transformation(origin = {52, 58},extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Modelica.Fluid.Sources.MassFlowSource_T Cold1(
      redeclare package Medium = Medium,
      T= 300,
      m_flow= 0.2,nPorts= 1, use_T_in = false, use_m_flow_in = false) annotation (Placement(visible = true, transformation(extent = {{-92, -10}, {-72, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp Ramp1(
      startTime=50,
      duration=5,
      height=0.4,
      offset=-0.2) annotation (Placement(visible = true, transformation(extent = {{76, 74}, {96, 94}}, rotation = 0)));
    inner Modelica.Fluid.System system(T_ambient (displayUnit = "K") = 303.15,energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyStateInitial,
        use_eps_Re=true) annotation (Placement(visible = true, transformation(extent = {{-94, 72}, {-74, 92}}, rotation = 0)));
  equation
//"constant";
    a = 10.98;
    b=-1.5472;
    ga=0.96*(10^(-10));
    E=22618;
    R=8.314472;
//"operation data"
  //Polynomials obtained by curve fitting with the operation data
    T1in = 439.913892 + 0.809553676 * (time / 86400) - 0.00438174426 * (time / 86400) ^ 2 + 0.0000116012692 * (time / 86400) ^ 3 - 1.59364124E-08 * (time / 86400) ^ 4 + 1.08121833E-11 * (time / 86400) ^ 5 - 2.86049304E-15 * (time / 86400) ^ 6;
    T1out=473.451134 + 0.639850142*(time/86400) - 0.00385198427*(time/86400)^2 + 0.0000109716761*(time/86400)^3 - 1.59140593E-08*(time/86400)^4 + 1.12607396E-11*(time/86400)^5 - 3.08179644E-15*(time/86400)^6;
    T2in=600.358404 - 0.107927146*(time/86400) + 0.00127465585*(time/86400)^2 - 0.00000469679366*(time/86400)^3 + 7.78173840E-09*(time/86400)^4 - 5.97561534E-12*(time/86400)^5 + 1.72703768E-15*(time/86400)^6;
    T2out=541.821838 - 0.219355978*(time/86400) + 0.00222018509*(time/86400)^2 - 0.00000767983863*(time/86400)^3 + 1.22561044E-08*(time/86400)^4 - 9.19522024E-12*(time/86400)^5 + 2.60872754E-15*(time/86400)^6;
    Qref=4563000;
//"datasheet"
    F = 0.9916;
    At=232.5493155;
    Ud=324.17;
    Re_t=34298;
    Ra=0.002979546;
    tao=30;
    Rf2=Ra*(1-exp(-(time/86400)/tao));
// "real equations"
    Delta1 = T2in - T1out;
    Delta2=(T2out-T1in);
    Delta=Delta1-Delta2;
    DeltaLmF=F*(Delta/(log(Delta1/Delta2)));
    Tf2=(((T2in+T2out)/2)-((T1in+T1out)/2));
    Tf=((T1in+T1out)/2)+0.55*Tf2;
    Uf=1/((1/Ud)+Rf2);
    Qf=Uf*At*DeltaLmF;
    Qloss=Qref-Qf;
    der(Eloss)=Qloss;
    VF=0.001879*Eloss/1000;
    e=Qf/Qref;
    der(Rf)=a*(Re_t^b)*exp(-E/(R*Tf))-ga*(Re_t^0.4) "This model does not fit with experimental data";
  connect(Ramp1.y, Hot2.m_flow_in) annotation(
      Line(points = {{97, 84}, {92.5, 84}, {92.5, 50}, {62, 50}}, color = {0, 0, 127}));
  connect(Hot2.ports[1], HEX.port_b2) annotation(
      Line(points = {{42, 58}, {26, 58}, {26, 28}, {28, 28}}, color = {255, 0, 0}, thickness = 0.5));
  connect(HEX.port_a2, End2.ports[1]) annotation(
      Line(points = {{-26, -26}, {-26, -60}, {-40, -60}}, color = {255, 0, 0}, thickness = 0.5));
  connect(Cold1.ports[1], HEX.port_a1) annotation(
      Line(points = {{-72, 0}, {-32, 0}, {-32, 0}, {-32, 0}}, color = {0, 0, 127}, thickness = 0.5));
  connect(HEX.port_b1, End1.ports[1]) annotation(
      Line(points = {{34, 0}, {76, 0}, {76, 0}, {76, 0}}, color = {0, 0, 127}, thickness = 0.5));
    annotation (experiment(StopTime=500, Tolerance=
            1e-005),
      Documentation(info="<html>
  <p>The simulation start in steady state with counterflow operation. At time t = 50, the mass flow rate on the secondary circuit is changed to a negative value in 5 seconds. After a transient, the heat exchanger operates in co-current flow.</p>
  <p><img src=\"modelica://Modelica/Resources/Images/Fluid/Examples/HeatExchanger.png\" alt=\"HeatExchanger.png\"/></p>
  </html>"));
  end HeatExchangerFoulingSimulation;

  package BaseClasses "Additional models for heat exchangers"
    extends Modelica.Icons.BasesPackage;

    model BasicHX "Simple heat exchanger model"
      outer Modelica.Fluid.System system "System properties";
      //General
      parameter Modelica.SIunits.Length length(min=0) "Length of flow path for both fluids";
      parameter Integer nNodes(min=1) = 2 "Spatial segmentation";
      parameter Modelica.Fluid.Types.ModelStructure modelStructure_1=Modelica.Fluid.Types.ModelStructure.av_vb
        "Determines whether flow or volume models are present at the ports"
        annotation(Evaluate=true, Dialog(tab="General",group="Fluid 1"));
      parameter Modelica.Fluid.Types.ModelStructure modelStructure_2=Modelica.Fluid.Types.ModelStructure.av_vb
        "Determines whether flow or volume models are present at the ports"
        annotation(Evaluate=true, Dialog(tab="General",group="Fluid 2"));
      replaceable package Medium_1 = Modelica.Media.Water.StandardWater constrainedby
        Modelica.Media.Interfaces.PartialMedium "Fluid 1"
                                                        annotation(choicesAllMatching, Dialog(tab="General",group="Fluid 1"));
      replaceable package Medium_2 = Modelica.Media.Water.StandardWater constrainedby
        Modelica.Media.Interfaces.PartialMedium "Fluid 2"
                                                        annotation(choicesAllMatching,Dialog(tab="General", group="Fluid 2"));
      parameter Modelica.SIunits.Area crossArea_1 "Cross sectional area" annotation(Dialog(tab="General",group="Fluid 1"));
      parameter Modelica.SIunits.Area crossArea_2 "Cross sectional area" annotation(Dialog(tab="General",group="Fluid 2"));
      parameter Modelica.SIunits.Length perimeter_1 "Flow channel perimeter" annotation(Dialog(tab="General",group="Fluid 1"));
      parameter Modelica.SIunits.Length perimeter_2 "Flow channel perimeter" annotation(Dialog(tab="General",group="Fluid 2"));
      final parameter Boolean use_HeatTransfer = true
        "= true to use the HeatTransfer_1/_2 model";
    // Heat transfer
      replaceable model HeatTransfer_1 =
          Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.IdealFlowHeatTransfer
        constrainedby
        Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.PartialFlowHeatTransfer
        "Heat transfer model" annotation(choicesAllMatching, Dialog(tab="General", group="Fluid 1", enable=use_HeatTransfer));
    
      replaceable model HeatTransfer_2 =
          Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.IdealFlowHeatTransfer
        constrainedby
        Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.PartialFlowHeatTransfer
        "Heat transfer model" annotation(choicesAllMatching, Dialog(tab="General", group="Fluid 2", enable=use_HeatTransfer));
    
      parameter Modelica.SIunits.Area area_h_1 "Heat transfer area" annotation(Dialog(tab="General",group="Fluid 1"));
      parameter Modelica.SIunits.Area area_h_2 "Heat transfer area" annotation(Dialog(tab="General",group="Fluid 2"));
     //Wall
      parameter Modelica.SIunits.Length s_wall(min=0) "Wall thickness"
        annotation (Dialog(group="Wall properties"));
      parameter Modelica.SIunits.ThermalConductivity k_wall
        "Thermal conductivity of wall material"
        annotation (Dialog(group="Wall properties"));
      parameter Modelica.SIunits.SpecificHeatCapacity c_wall
        "Specific heat capacity of wall material"
        annotation(Dialog(tab="General", group="Wall properties"));
      parameter Modelica.SIunits.Density rho_wall "Density of wall material"
        annotation(Dialog(tab="General", group="Wall properties"));
      final parameter Modelica.SIunits.Area area_h=(area_h_1 + area_h_2)/2
        "Heat transfer area";
      final parameter Modelica.SIunits.Mass m_wall=rho_wall*area_h*s_wall "Wall mass";
    // Assumptions
      parameter Boolean allowFlowReversal = system.allowFlowReversal
        "allow flow reversal, false restricts to design direction (port_a -> port_b)"
        annotation(Dialog(tab="Assumptions"), Evaluate=true);
      parameter Modelica.Fluid.Types.Dynamics energyDynamics=system.energyDynamics
        "Formulation of energy balance"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
      parameter Modelica.Fluid.Types.Dynamics massDynamics=system.massDynamics
        "Formulation of mass balance"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
      parameter Modelica.Fluid.Types.Dynamics momentumDynamics=system.momentumDynamics
        "Formulation of momentum balance, if pressureLoss options available"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
    //Initialization pipe 1
      parameter Modelica.SIunits.Temperature Twall_start "Start value of wall temperature"
                                                                            annotation(Dialog(tab="Initialization", group="Wall"));
      parameter Modelica.SIunits.Temperature dT "Start value for tube_side.T - shell_side.T"
        annotation (Dialog(tab="Initialization", group="Wall"));
      parameter Boolean use_T_start=true
        "Use T_start if true, otherwise h_start"
        annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium_1.AbsolutePressure p_a_start1=Medium_1.p_default
        "Start value of pressure"
        annotation(Dialog(tab = "Initialization", group = "Fluid 1"));
      parameter Medium_1.AbsolutePressure p_b_start1=Medium_1.p_default
        "Start value of pressure"
        annotation(Dialog(tab = "Initialization", group = "Fluid 1"));
      parameter Medium_1.Temperature T_start_1=if use_T_start then Medium_1.
          T_default else Medium_1.temperature_phX(
            (p_a_start1 + p_b_start1)/2,
            h_start_1,
            X_start_1) "Start value of temperature"
        annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 1", enable = use_T_start));
      parameter Medium_1.SpecificEnthalpy h_start_1=if use_T_start then Medium_1.specificEnthalpy_pTX(
            (p_a_start1 + p_b_start1)/2,
            T_start_1,
            X_start_1) else Medium_1.h_default
        "Start value of specific enthalpy"
        annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 1", enable = not use_T_start));
      parameter Medium_1.MassFraction X_start_1[Medium_1.nX]=Medium_1.X_default
        "Start value of mass fractions m_i/m"
        annotation (Dialog(tab="Initialization", group = "Fluid 1", enable=(Medium_1.nXi > 0)));
      parameter Medium_1.MassFlowRate m_flow_start_1 = system.m_flow_start
        "Start value of mass flow rate" annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 1"));
      //Initialization pipe 2
      parameter Medium_2.AbsolutePressure p_a_start2=Medium_2.p_default
        "Start value of pressure"
        annotation(Dialog(tab = "Initialization", group = "Fluid 2"));
      parameter Medium_2.AbsolutePressure p_b_start2=Medium_2.p_default
        "Start value of pressure"
        annotation(Dialog(tab = "Initialization", group = "Fluid 2"));
      parameter Medium_2.Temperature T_start_2=if use_T_start then Medium_2.
          T_default else Medium_2.temperature_phX(
            (p_a_start2 + p_b_start2)/2,
            h_start_2,
            X_start_2) "Start value of temperature"
        annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 2", enable = use_T_start));
      parameter Medium_2.SpecificEnthalpy h_start_2=if use_T_start then Medium_2.specificEnthalpy_pTX(
            (p_a_start2 + p_b_start2)/2,
            T_start_2,
            X_start_2) else Medium_2.h_default
        "Start value of specific enthalpy"
        annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 2", enable = not use_T_start));
      parameter Medium_2.MassFraction X_start_2[Medium_2.nX]=Medium_2.X_default
        "Start value of mass fractions m_i/m"
        annotation (Dialog(tab="Initialization", group = "Fluid 2", enable=Medium_2.nXi>0));
      parameter Medium_2.MassFlowRate m_flow_start_2 = system.m_flow_start
        "Start value of mass flow rate" annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 2"));
    //Pressure drop and heat transfer
      replaceable model FlowModel_1 =
          Modelica.Fluid.Pipes.BaseClasses.FlowModels.DetailedPipeFlow
        constrainedby
        Modelica.Fluid.Pipes.BaseClasses.FlowModels.PartialStaggeredFlowModel
        "Characteristic of wall friction" annotation(choicesAllMatching, Dialog(tab="General", group="Fluid 1"));
      replaceable model FlowModel_2 =
          Modelica.Fluid.Pipes.BaseClasses.FlowModels.DetailedPipeFlow
        constrainedby
        Modelica.Fluid.Pipes.BaseClasses.FlowModels.PartialStaggeredFlowModel
        "Characteristic of wall friction" annotation(choicesAllMatching, Dialog(tab="General", group="Fluid 2"));
      parameter Modelica.Fluid.Types.Roughness roughness_1=2.5e-5
        "Absolute roughness of pipe (default = smooth steel pipe)" annotation(Dialog(tab="General", group="Fluid 1"));
      parameter Modelica.Fluid.Types.Roughness roughness_2=2.5e-5
        "Absolute roughness of pipe (default = smooth steel pipe)" annotation(Dialog(tab="General", group="Fluid 2"));
    //Display variables
      Modelica.SIunits.HeatFlowRate Q_flow_1 "Total heat flow rate of pipe 1";
      Modelica.SIunits.HeatFlowRate Q_flow_2 "Total heat flow rate of pipe 2";
    
    HeatExchangerFouling5.BaseClasses.WallConstProps wall(
        rho_wall=rho_wall,
        c_wall=c_wall,
        T_start=Twall_start,
        k_wall=k_wall,
        dT=dT,
        s=s_wall,
        energyDynamics=energyDynamics,
        n=nNodes,
        area_h=area_h)
        annotation (Placement(transformation(extent={{-29,-23},{9,35}})));
    
    Modelica.Fluid.Pipes.DynamicPipe tube_side(
        redeclare final package Medium = Medium_1,
        final isCircular=false,
        final diameter=0,
        final nNodes=nNodes,
        final allowFlowReversal=allowFlowReversal,
        final energyDynamics=energyDynamics,
        final massDynamics=massDynamics,
        final momentumDynamics=momentumDynamics,
        final length=length,
        final use_HeatTransfer=use_HeatTransfer,
        redeclare final model HeatTransfer = HeatTransfer_1,
        final use_T_start=use_T_start,
        final T_start=T_start_1,
        final h_start=h_start_1,
        final X_start=X_start_1,
        final m_flow_start=m_flow_start_1,
        final perimeter=perimeter_1,
        final crossArea=crossArea_1,
        final p_a_start=p_a_start1,
        final p_b_start=p_b_start1,
        final roughness=roughness_1,
        redeclare final model FlowModel = FlowModel_1,
        final modelStructure=modelStructure_1) annotation (Placement(visible = true, transformation(extent = {{-40, -70}, {20, -10}}, rotation = 0)));
    
    Modelica.Fluid.Pipes.DynamicPipe shell_side(
        redeclare final package Medium = Medium_2,
        final nNodes=nNodes,
        final allowFlowReversal=allowFlowReversal,
        final energyDynamics=energyDynamics,
        final massDynamics=massDynamics,
        final momentumDynamics=momentumDynamics,
        final length=length,
        final isCircular=false,
        final diameter=0,
        final use_HeatTransfer=use_HeatTransfer,
        redeclare final model HeatTransfer = HeatTransfer_2,
        final use_T_start=use_T_start,
        final T_start=T_start_2,
        final h_start=h_start_2,
        final X_start=X_start_2,
        final m_flow_start=m_flow_start_2,
        final perimeter=perimeter_2,
        final crossArea=crossArea_2,
        final p_a_start=p_a_start2,
        final p_b_start=p_b_start2,
        final roughness=roughness_2,
        redeclare final model FlowModel = FlowModel_2,
        final modelStructure=modelStructure_2)
                  annotation (Placement(visible = true, transformation(extent = {{20, 80}, {-40, 20}}, rotation = 0)));
    
      Modelica.Fluid.Interfaces.FluidPort_b port_b1(redeclare final package
          Medium =
            Medium_1) annotation (Placement(transformation(extent={{100,-12},{120,
                8}})));
      Modelica.Fluid.Interfaces.FluidPort_a port_a1(redeclare final package
          Medium =
            Medium_1) annotation (Placement(transformation(extent={{-120,-12},{
                -100,8}})));
      Modelica.Fluid.Interfaces.FluidPort_b port_b2(redeclare final package
          Medium =
            Medium_2) annotation (Placement(visible = true, transformation(extent = {{80, 86}, {100, 106}}, rotation = 0), iconTransformation(extent = {{80, 80}, {100, 100}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_a port_a2(redeclare final package
          Medium =
            Medium_2) annotation (Placement(visible = true, transformation(extent = {{-100, -102}, {-80, -82}}, rotation = 0), iconTransformation(extent = {{-100, -100}, {-80, -80}}, rotation = 0)));
    
    equation
      Q_flow_1 = sum(tube_side.heatTransfer.Q_flows);
      Q_flow_2 = sum(shell_side.heatTransfer.Q_flows);
    connect(shell_side.port_b, port_b2) annotation(
        Line(points = {{-40, 50}, {-76, 50}, {-76, 96}, {90, 96}}, color = {255, 0, 0}, thickness = 0.5));
    connect(tube_side.port_b, port_b1) annotation(
        Line(points = {{20, -40}, {42, -40}, {42, -2}, {110, -2}}, color = {0, 0, 127}, thickness = 0.5));
    connect(tube_side.port_a, port_a1) annotation(
        Line(points = {{-40, -40}, {-75.3, -40}, {-75.3, -2}, {-110, -2}}, color = {0, 0,127}, thickness = 0.5));
    connect(shell_side.port_a, port_a2) annotation(
        Line(points = {{20, 50}, {65, 50}, {65, -92}, {-90, -92}}, color = {255, 0, 0}, thickness = 0.5));
    connect(wall.heatPort_b, tube_side.heatPorts) annotation(
        Line(points = {{-10, -8.5}, {-10, -27}}, color = {191, 0, 0}));
    connect(shell_side.heatPorts[nNodes:(-1):1], wall.heatPort_a[1:nNodes]) annotation(
        Line(points = {{-10, 37}, {-10, 20.5}}, color = {127, 0, 0}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, initialScale = 0.1), graphics={Rectangle(fillColor = {95, 95, 95}, fillPattern = FillPattern.Forward, extent = {{-100, -26}, {100, -30}}), Rectangle(fillColor = {95, 95, 95}, fillPattern = FillPattern.Forward, extent = {{-100, 30}, {100, 26}}), Rectangle(fillColor = {0, 63, 125}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, 60}, {100, 30}}), Rectangle(fillColor = {0, 63, 125}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, -30}, {100, -60}}), Rectangle(fillColor = {0, 128, 255}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, 26}, {100, -26}}), Text(origin = {-2, -86}, lineColor = {255, 255, 255}, extent = {{-150, 110}, {150, 70}}, textString = "%name"), Polygon(origin = {-20, 16},lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid, points = {{20, -80}, {40, -89}, {20, -100}, {20, -80}}), Polygon(origin = {50, 0},lineColor = {255, 0, 0}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, points = {{-50, 82}, {-70, 73}, {-50, 62}, {-50, 82}}), Rectangle(origin = {90, 70}, fillColor = {0, 63, 125}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-10, 10}, {10, -10}}), Rectangle(origin = {-90, -70}, fillColor = {0, 63, 125}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-10, 10}, {10, -10}})}),
        Documentation(info="<html>
    <p>Simple model of a heat exchanger consisting of two pipes and one wall in between.
    For both fluids geometry parameters, such as heat transfer area and cross section as well as heat transfer and pressure drop correlations may be chosen.
    The flow scheme may be concurrent or counterflow, defined by the respective flow directions of the fluids entering the component.
    The design flow direction with positive m_flow variables is counterflow.</p>
    </html>"));
    end BasicHX;

    model WallConstProps "Pipe wall with capacitance, assuming 1D heat conduction and constant 
    material properties"
    Modelica.SIunits.ThermalInsulance Rf "Fouling resistance";
      parameter Integer n(min = 1) = 1 "Segmentation perpendicular to heat conduction";
      //Geometry
      parameter Modelica.SIunits.Length s "Wall thickness";
      parameter Modelica.SIunits.Area area_h "Heat transfer area";
      //Material properties
      parameter Modelica.SIunits.Density rho_wall "Density of wall material";
      parameter Modelica.SIunits.SpecificHeatCapacity c_wall "Specific heat capacity of wall material";
      parameter Modelica.SIunits.ThermalConductivity k_wall "Thermal conductivity of wall material";
      parameter Modelica.SIunits.Mass[n] m = fill(rho_wall * area_h * s / n, n) "Distribution of wall mass";
      //Initialization
      outer Modelica.Fluid.System system;
      parameter Modelica.Fluid.Types.Dynamics energyDynamics = system.energyDynamics "Formulation of energy balance" annotation(
        Evaluate = true,
        Dialog(tab = "Assumptions", group = "Dynamics"));
      parameter Modelica.SIunits.Temperature T_start "Wall temperature start value";
      parameter Modelica.SIunits.Temperature dT "Start value for port_b.T - port_a.T";
      //Temperatures
      Modelica.SIunits.Temperature[n] Tb(each start = T_start + 0.5 * dT);
      Modelica.SIunits.Temperature[n] Ta(each start = T_start - 0.5 * dT);
      Modelica.SIunits.Temperature[n] T(start = ones(n) * T_start, each stateSelect = StateSelect.prefer) "Wall temperature";
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] heatPort_a "Thermal port" annotation(
        Placement(visible = true, transformation(extent = {{-20, 40}, {20, 60}}, rotation = 0), iconTransformation(extent = {{-20, 60}, {20, 80}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] heatPort_b "Thermal port" annotation(
        Placement(visible = true, transformation(extent = {{-20, -40}, {20, -60}}, rotation = 0), iconTransformation(extent = {{-20, -60}, {20, -80}}, rotation = 0)));
    initial equation
      if energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyStateInitial then
        der(T) = zeros(n);
      elseif energyDynamics == Modelica.Fluid.Types.Dynamics.FixedInitial then
        T = ones(n) * T_start;
      end if;
    equation
    der(Rf)=1000*exp(-time);
      for i in 1:n loop
        assert(m[i] > 0, "Wall has negative dimensions");
        if energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState then
          0 = heatPort_a[i].Q_flow + heatPort_b[i].Q_flow;
        else
          c_wall * m[i] * der(T[i]) = heatPort_a[i].Q_flow + heatPort_b[i].Q_flow;
        end if;
        heatPort_a[i].Q_flow = 2 * k_wall / s * (Ta[i] - T[i]) * area_h / n;
        heatPort_b[i].Q_flow = 2 * k_wall / s * (Tb[i] - T[i]) * area_h / n;
      end for;
      Ta = heatPort_a.T;
      Tb = heatPort_b.T;
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1), graphics = {Rectangle(origin = {0, 20}, fillColor = {95, 95, 95}, fillPattern = FillPattern.Forward,extent = {{-100, 40}, {100, -40}}), Text(origin = {-2, 20},extent = {{-82, 18}, {76, -18}}, textString = "%name"), Rectangle(origin = {0, -40}, lineColor = {140, 70, 0}, fillColor = {223, 223, 223}, fillPattern = FillPattern.CrossDiag, extent = {{-100, 20}, {100, -20}}), Text(origin = {2, -36}, extent = {{-82, 18}, {76, -18}}, textString = "Fouling")}),
        Documentation(revisions = "<html>
    <ul>
    <li><em>04 Mar 2006</em>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
    </ul>
    </html>", info = "<html>
    Simple model of circular (or any other closed shape) wall to be used for pipe (or duct) models. Heat conduction is regarded one dimensional, capacitance is lumped at the arithmetic mean temperature. The spatial discretization (parameter <code>n</code>) is meant to correspond to a connected fluid model discretization.
    </html>"));
    end WallConstProps;
  end BaseClasses;
end HeatExchangerFouling5;
