# calc_pulse_tool
This tool uses  **[Brian Hargreave's Bloch Simulator](http://mrsrl.stanford.edu/~brian/bloch/)** to perform Bloch simulations of the RF Pulses in Matlab. Also has the option to simulate inhomogeneous B1 excitation fields as well as repeated excitations (bSSFP).

## Installation:
1. **[Brian Hargreave's Bloch Simulator](http://mrsrl.stanford.edu/~brian/bloch/)** 
2. You'll need to have Matlab's Appdesigner installed.

## Usage:
1. Start by either double-clicking _simulate_rf_pulse_b1_tool.mlapp_ (should start appdesigner) or by typing appdesigner in Matlab's command window and opening the _simulate_rf_pulse_b1_tool.mlapp_
2. Run the tool by hitting function key F5 or pressing the green "Play" button.
3. Load an RF Pulse (either .mat or .dat in the _RF Pulses_ subfolder)
    1. slr_sharp1_1ms_2000Hz will give you a Shinnar-LeRouge Pulse of 1ms and 2000Hz resulting bandwidth with "sharpness 1" (taken from Bruker).
    2. AFP_2_3182ms_Hsn is an Adibatic Full Passage Pulse of 2.3162ms (taken from Robin de Graaf's pulsewizard toolbox)
4. Press load to load the pulse.
5. You can now edit pulse duration, flipangle, pulse amplitude, and gradient strength. Bandwith factor BWfac and integration factor Sint are intrisic properties of the pulse and can not be changed.
You can also change the range and the resolution (frequency and position) over which the simulation will be performed. 
6. Pressing *Simulate!* will run the Bloch Simulation with the selected parameters.
7. The simulation results will be plotted in the lower left windows.
1. You can select one or multiple (hold Shift key) properties (Mx, My, ...) of the simulated magnetization.
2. You can select multiple frequencies by choosing the *Show 2 Freqs*. Then you can select Freq 1 and Freq 2.
3. You can also use the slider to select different frequencies.
8. On the Bottom you can find a filed saying "save". For now I'd recommend to select all options on the list and hit the _save_ button.



