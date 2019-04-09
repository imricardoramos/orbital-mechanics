try
    componentAstrogator = scenario.ComponentDirectory.GetComponents('eComponentAstrogator');
catch
    app = actxserver('STK11.application');
    app.UserControl = 1;
    app.visible = 1;
    root = app.Personality2;
    root.NewScenario('Validation');
    scenario = root.CurrentScenario;

    scenario.SetTimePeriod('24 Feb 2012 16:00:00.000','17 Sep 2012 17:00:00.000');
    scenario.StartTime = '24 Feb 2012 16:00:00.000';
    scenario.StopTime = '17 Sep 2012 17:00:00.000';

    root.ExecuteCommand('Animate * Reset');
    %% Satellite Definition and Astrogator Init
    try
        scenario.Children.Unload('eSatellite','ValidationSat');
    catch
    end
    satellite = scenario.Children.New('eSatellite','ValidationSat');
    satellite.SetPropagatorType('ePropagatorAstrogator');
    ASTG = satellite.Propagator;
    MCS = ASTG.MainSequence;
end
%% Astrogator Config
componentAstrogator = scenario.ComponentDirectory.GetComponents('eComponentAstrogator');
engineModels = componentAstrogator.GetFolder('Engine Models');
try
    customThruster = engineModels.Item('customThruster');
catch
    customThruster = engineModels.DuplicateComponent('Constant Thrust and Isp','customThruster');
end
customThruster.Isp = 720;
customThruster.Thrust = 0.001;

MCS.RemoveAll;

initstate = MCS.Insert('eVASegmentTypeInitialState','Inner Orbit','-');
initstate.OrbitEpoch = scenario.StartTime;
initstate.SetElementType('eVAElementTypeModKeplerian');

initstate.InitialState.DryMass = 0.10000000000000009;
initstate.InitialState.Cd = 2.2;
initstate.InitialState.DragArea = 10;
initstate.InitialState.Cr = 1;
initstate.InitialState.SRPArea = 20;

initstate.FuelTank.FuelMass = 2.9;

initstate.Element.RadiusOfPeriapsis = 7000.0;
initstate.Element.Eccentricity = 0.06666666666666667;
initstate.Element.Inclination = 65.1;
initstate.Element.RAAN = 340.0;
initstate.Element.ArgOfPeriapsis = 58.0;
initstate.Element.TrueAnomaly = 332.00000000000006;

propagate = MCS.Insert('eVASegmentTypePropagate','Propagate','-');
propagate.PropagatorName = 'Earth Point Mass';
propagate.Properties.Color = uint32(hex2dec('00ff00'));
propagate.StoppingConditions.Item('Duration').Properties.Trip = 3600;

dv1 = MCS.Insert('eVASegmentTypeManeuver','DV1','-');
dv1.Properties.Color = uint32(hex2dec('00d3ff'));
dv1.SetManeuverType('eVAManeuverTypeFinite');
maneuver = dv1.Maneuver;
maneuver.SetAttitudeControlType('eVAAttitudeControlThrustVector');
maneuver.AttitudeControl.ThrustAxesName = 'Satellite VNC(Earth)';
maneuver.AttitudeControl.ThrustVector.AssignXYZ(1,0,0);
maneuver.Propagator.StoppingConditions.Item('Duration').Properties.Trip = 6307200.0;
maneuver.Propagator.PropagatorName = 'Earth Point Mass';
maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','customThruster');

dv2 = MCS.Insert('eVASegmentTypeManeuver','DV2','-');
dv2.Properties.Color = uint32(hex2dec('00ffd3'));
dv2.SetManeuverType('eVAManeuverTypeFinite');
maneuver2 = dv2.Maneuver;
maneuver2.SetAttitudeControlType('eVAAttitudeControlThrustVector');
maneuver2.AttitudeControl.ThrustAxesName = 'Satellite VNC(Earth)';
maneuver2.AttitudeControl.ThrustVector.AssignXYZ(0,1,0);
maneuver2.Propagator.StoppingConditions.Item('Duration').Properties.Trip = 6307200.0;
maneuver2.Propagator.PropagatorName = 'Earth Point Mass';
maneuver2.SetPropulsionMethod('eVAPropulsionMethodEngineModel','customThruster');

propagate2 = MCS.Insert('eVASegmentTypePropagate','Propagate','-');
propagate2.PropagatorName = 'Earth Point Mass';
propagate2.Properties.Color = uint32(hex2dec('00ff00'));
propagate2.StoppingConditions.Item('Duration').Properties.Trip = 5184000;

ASTG.RunMCS;

keplerianElemsDP = satellite.DataProviders.Item('Astrogator Values').Group.Item('Keplerian Elems').Exec(scenario.StartTime, scenario.StopTime,60);
maneuverDP = satellite.DataProviders.Item('Astrogator Values').Group.Item('Maneuver').Exec(scenario.StartTime, scenario.StopTime,60);
keplerianElemsData = keplerianElemsDP.DataSets.ToArray;
maneuverData = maneuverDP.DataSets.ToArray;
maneuverData = maneuverData(:,2:end);
table = cell2table([keplerianElemsData,maneuverData],'VariableNames',[keplerianElemsDP.DataSets.ElementNames;maneuverDP.DataSets.ElementNames(2:end)]);
writetable(table,'STKValidationScriptMatlab_DataOut.csv');

exit;