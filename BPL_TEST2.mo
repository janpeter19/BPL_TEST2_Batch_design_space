encapsulated package BPL_TEST2

	// Test culture with 2 states - print font size: 8 
	//
	// Author: Jan Peter Axelsson
	// 2019-02-02 - Created
	// 2019-02-17 - Renewed work on dividing code in library and application
	// 2019-02-17 - Adapt also EnvironmentCon similar to LiquidCon
	// 2019-02-18 - Change Medium.X=1 from 2 as a built in constant
	// 2019-03-11 - Use LiquidCon and EnvironmentCon as parameters for package Equipment
	// 2019-03-12 - Replace EnvironmentCon with inner/outer connection
	// 2019-03-13 - Works!
	// 2019-03-14 - Corrected import of package Equipment - works
	// 2019-05-21 - Many changes following DEMO_v22 style and introduced record Medium2_data
	// 2019-05-22 - In ExploreCultureType I address MediumCHO.nc now, since not imported
	// 2019-05-24 - Try to make configurations general code and change style of names with that purpose in mind
	// 2020-02-19 - Works for both JModelica 2.14 and OpenModelica 1.16 nightly build
	// 2020-03-16 - Adapted for BR5m with stream and pressure etc
	// 2020-07-13 - Start as BP6a_TEST2 and imported from BR5m
	// 2020-07-29 - Changed handling of nc and now what I call functional approach of adaptation
	// 2020-10-01 - Adapt to development of BP6a from BP6c
	// 2020-10-10 - Simplified Yxs to Y 
	// 2020-11-05 - Introduced variable D for dilution or perfusion rate for Chemostat and Perfusion setup
	// 2020-11-21 - Adapt for ReactorType with arbitrary number of inlets, outlets and ports
	// 2020-11-25 - Final adaptation for general ReactorTypes
	// 2021-02-03 - Changed MediumBase to LiquidphaseBase etc
	// 2021-02-08 - Changed where appropriate to declaration Liquidphase.Concentration c
	// 2021-02-10 - Change file name to BPL and drop version number - towards BPL_v2.0
	// 2021-02-12 - Introduced package Interfaces to cover media and ReactorInterface
	// 2021-02-13 - Introduced package Media to cover at the moment only standard air
	// 2021-02-16 - Introduced package UsersGuide to bring BPL version etc
	// 2021-02-19 - Included import BPL.Interfaces... in application package... encapsulated packge require that
	// 2021-02-20 - Changed to media types also for mass and rate
	// 2021-02-21 - Start to eliminate un-necessary suffix Type
	// 2021-02-23 - Finished elimination of un-necessary suffix Type
	// 2021-02-23 - Introduced new adaptation of EquipmentLib and Reactor - idea from Adrian Pop
	// 2021-02-26 - I corrected model ExploreCultureSystem concerning adaptation to Liquidphase2
	// 2021-03-18 - I corrected the modification of the library - Stephane at Modelon saw it
	// 2021-03-19 - Change of inner/outer connection - Stéphande Veluts recomendation
	// 2021-03-20 - Corrected use of import from ReactorInterfaceInner to Reactor and adaptation EquipmentLib
	// 2021-04-22 - Tidy up indentations and introduced NoVolumeRate in adaptation of EquipmentLib
	// 2021-05-02 - Simplified adaptation of Reactor after small change in BPL 2.0.5 
	// 2021-07-15 - Included a setup for batch expansion from reactor bio1 to bio2 both in batch mode
	// 2021-11-22 - Corrected import from MSL for Sources and Types and now pass OpenModelica 1.19.0
	// 2022-05-28 - Introduced a variable mu in the culture to simplify access and improve transparency for user
	// 2022-08-18 - Introduce MSL_info and include in all models
	// 2022-08-19 - Moved MSL_info to BPL.UsersGuide
	// 2022-09-03 - New file TEST2B but still call the package BPL_TEST2 and extend with measurement noise etc
	// 2022-09-04 - Refined declaration of MSL.usage to reflect what components are used for each sytem
	// 2022-09-06 - I see the need for a separate system without noise to avoid sampling and rely on events
 
// ----------------------------------------------------------------------------------------------------------------
//    Customer informtion
// ---------------------------------------------------------------------------------------------------------------- 
   
	record Customer_info
		import BPL.UsersGuide.Customer_info;
		extends Customer_info(
		   name = "Bioprocess Library", 
		   address = "Drejargatan 1, 113 42 Stockholm, Sweden",
			contact_person = "Jan Peter Axelsson", 
			license_valid_until = {2022, 12, 31});	
	end Customer_info;

// ----------------------------------------------------------------------------------------------------------------
//    Define Liquidphase2 and Culture2
// ----------------------------------------------------------------------------------------------------------------
   
   package Liquidphase2
		import BPL.Interfaces.LiquidphaseBase;
		extends LiquidphaseBase(
			name="Standard components X and S",
			nc=2);		
		constant Integer X=1                                     "Index cell";
		constant Integer S=2                                     "Index substrate"; 
		constant Real[nc] mw (each unit="Da")  
		                   = {24.6, 180.0}                       "Molecular weight";  
	end Liquidphase2;

	record Liquidphase_data
		constant String name = Liquidphase2.name;
		constant Integer nc = Liquidphase2.nc;
		constant Integer X  = Liquidphase2.X                     "Cells";
		constant Integer S = Liquidphase2.S                      "Substrate";
		constant Real[nc] mw = Liquidphase2.mw;
	end Liquidphase_data;

	model Culture2
		// Connection of reactor micro-environment concentrations c to specific cell flows q:
		import BPL.Interfaces.ReactorInterface;
		extends ReactorInterface;
		outer Liquidphase.Concentration c;
		Liquidphase.Rate q;	

		// Liquidphase constants:
		import BPL_TEST2.Liquidphase2.*;
	
		// Parameters
      parameter Real Ks (unit="g/L") = 0.1                     "Substrate uptake saturation";
      parameter Real qSmax (unit="g/(g*h)") = 1.00             "Substrate uptake maximal rate";
      parameter Real Y (unit="g/g") = 0.50                     "Yield of cells from substrate";
		// Variable
		Real mu (unit="1/h")                                     "Cell specific growth rate variable";       	
	equation
		q[S] = - qSmax*c[S]/(c[S]+Ks);
		mu = -Y*q[S];
		q[X] = mu;
	end Culture2;

// ----------------------------------------------------------------------------------------------------------------
//    Define a testbed to explore the culture response from c[] to flow q[]
// ----------------------------------------------------------------------------------------------------------------

	model ExploreCultureSystem
		import BPL.Interfaces.ReactorInterface;                                
		extends ReactorInterface(redeclare package Liquidphase=Liquidphase2);
		
		// Actual culture
		Culture2 culture(redeclare package Liquidphase=Liquidphase2);	
		
		// Connection to culture
		inner Liquidphase.Concentration c;
		
		// Constants data about the medium:
		import BPL_TEST2.Liquidphase2.*;
		
		// Table of concentrations
		Integer index;
		parameter Integer n = 50;
		parameter Real[n] table = 0.01:0.01:0.50;
				
		// Sampling
		parameter Real samplePeriod = 1;
		parameter Real sampleStart = 0;
	equation
		when sample(sampleStart, samplePeriod) then
			index = if time < 1 then 1 else integer(ceil(time));
			c[X] = 1;
			c[S] = table[index];			
		end when;
	end ExploreCultureSystem;

// ----------------------------------------------------------------------------------------------------------------
//    Adaptation of the library Equipment to the actual medium and culture and add two blocks
// ----------------------------------------------------------------------------------------------------------------

	import BPL.Interfaces.ReactorInterfaceInner;
 	import Modelica.Blocks.Noise.NormalNoise;

	package Equipment
		import BPL.EquipmentLib;
	   extends EquipmentLib(
			redeclare package Liquidphase = Liquidphase2,
			Reactor(redeclare model Culture = Culture2(redeclare package Liquidphase=Liquidphase2)));

		block SampleNoise
			LiquidCon probe;
			discrete output LiquidCon out;
			import BPL_TEST2.Liquidphase2.*;                        // Get mnemomnics X and S  
			parameter Real sigma (unit="g/L") = 0.1                 "Standard deviation on measured substrate conc S";
			parameter Real samplePeriod (unit="h") = 0.1            "Sample period of noise generator";
			parameter Real sampleStart (unit="h") = samplePeriod    "Start time";
			Real p (unit="bar")                                     "Pressure";	

			NormalNoise noise(samplePeriod=samplePeriod, mu=0.0, sigma=sigma);
				inner Modelica.Blocks.Noise.GlobalSeed globalSeed;   	
		equation
			probe.p = p;
			out.p = p;
			probe.F = 0;
//			out.F = 0;
			when sample(sampleStart, samplePeriod) then
				out.c[X] = probe.c[X];
				out.c[S] = probe.c[S] + noise.y;
			end when;
			for i in 1:Liquidphase.nc loop
				inStream(probe.c[i]) = probe.c[i];
			end for;
		end SampleNoise;
	   	
		block DetectEndBatch
			LiquidCon probe;
			import BPL_TEST2.Liquidphase2.*;                        // Get mnemomnics X and S         
			parameter Real S_min (unit="g/L") = 1.0                 "Substrate conc limitation def end of batch";
			parameter Real time_final_max (unit="h") = 6.0          "Specification of maximal time_final value";
			parameter Real X_final_min (unit="g/L") = 5.0           "Specification of minimal X_final value";
			Real time_final (start=0, fixed=true, unit="h")         "Time final, when substrate conc goes below S_min";
			Real X_final (start=0, fixed=true, unit="g/L")          "Cell conc final";
			Real S_final (start=0, fixed=true, unit="g/L")          "Substrate conc final";
			Real p (unit="bar")                                     "Pressure";
			Real batch_evaluation (start=0, fixed=true)             "Batch evaluation - accepted for >0";
		   discrete Boolean firstTime (start=true, fixed=true)     "Detect crossing S<S_min firstTime only";
		equation
			probe.p = p;
			probe.F = 0;
			when (probe.c[S] < S_min) and pre(firstTime) then
				firstTime = false;
				time_final = time;
				X_final = inStream(probe.c[X]);
				S_final = inStream(probe.c[S]);
				batch_evaluation = if (time < time_final_max) and (inStream(probe.c[X]) > X_final_min) then 1 
								   elseif (time > time_final_max) and (inStream(probe.c[X]) > X_final_min) then -1
								   elseif (time < time_final_max) and (inStream(probe.c[X]) < X_final_min) then -2
								   else -3;
			end when;
			for i in 1:Liquidphase.nc loop
				inStream(probe.c[i]) = probe.c[i];
			end for;
		end DetectEndBatch;					
	end Equipment;	
	
// ----------------------------------------------------------------------------------------------------------------
//    Connecting systems for different examples of operation
// ----------------------------------------------------------------------------------------------------------------
	
	import BPL_TEST2.Equipment;
	import BPL.UsersGuide.BPL_info;
	import BPL.UsersGuide.MSL_info;
 	import Modelica.Blocks.Sources; 
	import Modelica.Blocks.Types;

	model ExploreCulture
		BPL_info BPL; 
		MSL_info MSL(usage=" ");
		Customer_info Customer;	
		Liquidphase_data liquidphase;
		ExploreCultureSystem testbed;
	end ExploreCulture;	

   model Batch "Batch cultivation"
		BPL_info BPL; 
		MSL_info MSL(usage=" ");
		Customer_info Customer;	
		Liquidphase_data liquidphase;
      Equipment.Reactor bioreactor(X=liquidphase.X);
   equation
   end Batch;	

	import BPL.Control;

   model Fedbatch "Fedbatch cultivation"
		BPL_info BPL; 
		MSL_info MSL(usage="RealÍnput, RealOutput");
		Customer_info Customer;	
		Liquidphase_data liquidphase;
      Equipment.Reactor bioreactor(X=liquidphase.X, n_inlets=1);
      Equipment.FeedSystem feedtank;
      Control.DosageSchemeExp dosagescheme;
   equation
      connect(bioreactor.inlet[1],feedtank.outlet);
      connect(dosagescheme.F,feedtank.Fsp);
   end Fedbatch;

	model Batch_expansion "Batch expansion"
		BPL_info BPL; 
		MSL_info MSL(usage="RealÍnput, RealOutput, CombiTimeTable, Types");
		Customer_info Customer;	
		Liquidphase_data liquidphase;
		Equipment.Reactor bio1(X=liquidphase.X, n_outlets=1);
		Equipment.Reactor bio2(X=liquidphase.X, n_inlets=1);	
		Equipment.Pump pump;	
      Sources.CombiTimeTable schemePump(
		   smoothness=Types.Smoothness.ConstantSegments,
			extrapolation = Types.Extrapolation.HoldLastPoint,
		   table=[0,0; 1001,1; 1002,2; 1003,3; 1004,4; 1005,5; 1006,6; 1007,7; 1008,8; 1009,9]);
		// Help variables for specific growth rates
		Real bio1_culture_mu	(unit="1/g")  "Bio1 relevant specific growth rate";
		Real bio2_culture_mu	(unit="1/g")  "Bio2 relevant specific growth rate";
	equation
		bio1_culture_mu = if bio1.m[1] > 0.01 then bio1.culture.q[1] else 0;
		bio2_culture_mu = if bio2.m[1] > 0.01 then bio2.culture.q[1] else 0;		
		connect(bio1.outlet[1], pump.inlet);
		connect(pump.outlet, bio2.inlet[1]);
		connect(schemePump.y[1], pump.Fsp);
	end Batch_expansion;

	model Batch_harvest "Batch cultivation with harvest using ideal filtration"
		BPL_info BPL; 
		MSL_info MSL(usage="RealÍnput, RealOutput, CombiTimeTable, Types");
		Customer_info Customer;	
		Liquidphase_data liquidphase;
		Equipment.Reactor bioreactor(X=liquidphase.X, n_outlets=1);
		Equipment.FilterSystem filter;
		Equipment.HarvestTank harvesttank_cells, harvesttank_liquid;
		Equipment.Pump pump;	
      Sources.CombiTimeTable schemePump(
		   smoothness=Types.Smoothness.ConstantSegments,
			extrapolation=Types.Extrapolation.HoldLastPoint,
		   table=[0,0; 1001,1; 1002,2; 1003,3; 1004,4; 1005,5; 1006,6; 1007,7; 1008,8; 1009,9]);
	equation
		connect(bioreactor.outlet[1], pump.inlet);
		connect(pump.outlet, filter.inlet);
		connect(filter.retentate, harvesttank_cells.inlet);
		connect(filter.filtrate, harvesttank_liquid.inlet);
		connect(schemePump.y[1], pump.Fsp);
	end Batch_harvest;

   model Perfusion "Perfusion cultivation"
		BPL_info BPL;
		MSL_info MSL(usage="RealÍnput, RealOutput, CombiTimeTable, Types");
		Customer_info Customer;	
		Liquidphase_data liquidphase;
      Equipment.Reactor bioreactor(X=liquidphase.X, n_inlets=2, n_outlets=1);
		Equipment.FilterSystem filter;
      Equipment.FeedSystem feedtank;
      Equipment.HarvestTank harvesttank;
      Sources.CombiTimeTable schemePump1(
		   smoothness=Types.Smoothness.ConstantSegments,
			extrapolation = Types.Extrapolation.HoldLastPoint,
		   table=[0,0; 1001,1; 1002,2; 1003,3; 1004,4; 1005,5; 1006,6; 1007,7; 1008,8; 1009,9]);	
      Sources.CombiTimeTable schemePump2(
		   smoothness=Types.Smoothness.ConstantSegments,
			extrapolation = Types.Extrapolation.HoldLastPoint,
		   table=[0,0; 1001,1; 1002,2; 1003,3; 1004,4; 1005,5; 1006,6; 1007,7; 1008,8; 1009,9]);			
		// Help variable for visualization of data
		Real D (unit="1/h") "Perfusion  rate";		
   equation
		D = bioreactor.inlet[1].F/bioreactor.V;
      connect(bioreactor.inlet[1], feedtank.outlet);
      connect(bioreactor.outlet[1], filter.inlet);
		connect(bioreactor.inlet[2], filter.retentate);
		connect(filter.filtrate, harvesttank.inlet);
      connect(schemePump1.y[1], feedtank.Fsp);
      connect(schemePump2.y[1], filter.Fsp);
   end Perfusion;

	model BatchWithNoise "Batch with end detection"
		BPL_info BPL; 
		MSL_info MSL(usage="Noise.NormalNoise");
		Customer_info Customer;	
		Liquidphase_data liquidphase;
      Equipment.Reactor bioreactor(X=liquidphase.X, n_ports=1);
		Equipment.SampleNoise sensor;
		Equipment.DetectEndBatch monitor;
   equation
		connect(bioreactor.port[1], sensor.probe);	
		connect(sensor.out, monitor.probe);
	end BatchWithNoise;

	model BatchNoNoise "Batch with end detection"
		BPL_info BPL; 
		MSL_info MSL(usage=" ");
		Customer_info Customer;	
		Liquidphase_data liquidphase;
      Equipment.Reactor bioreactor(X=liquidphase.X, n_ports=1);
		Equipment.DetectEndBatch monitor;
   equation
		connect(bioreactor.port[1], monitor.probe);
	end BatchNoNoise;

end BPL_TEST2;