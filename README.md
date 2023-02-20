# calc_pulse_tool
This tool uses  **[Brian Hargreave's Bloch Simulator](http://mrsrl.stanford.edu/~brian/blochsim/)** to perform Bloch simulations of the RF Pulses in Matlab. Also has the option to simulate inhomogeneous B1 excitation fields as well as repeated excitations (bSSFP).

## Installation:
1. Install **[Brian Hargreave's Bloch Simulator](http://mrsrl.stanford.edu/~brian/blochsim/)**.
- Download *bloch.m* and *bloch.c* and save them into the folder *bloch_sim* (make a new one)
- Add the path where you store them to your matlab PATHS (HOME &rarr; Set Path &rarr; Add Folder)
- Open `bloch.c` with any editor and change 
1. all `malloc` to `mxMalloc`
2. all `free` to `mxFree`
3. in line **551** `int outsize[3];`  to `mwSize outsize[3];`

In case you don't want DEBUG output into your matlab Command Window, change `#define DEBUG` to `//#define DEBUG` line 8.


- Open a Matlab session in the *bloch_sim* directory or use `cd('\path\to\bloch_sim\')`
- In the Command Window run
    ```
    mex -v bloch.c
    ```
    This should result in
    ```
    MEX completed successfully.
    ```
- You might need to install an compiler (Windows), the one that MathWorks suggests is this one: https://de.mathworks.com/support/requirements/supported-compilers.html

2. Download `mri_rf_pulse_sim_tool` and add it to your matlabs PATHS (HOME &rarr; Set Path &rarr; Add Folder)
- You'll need to have Matlab's Appdesigner installed.
- type `appdesigner` in matlab Command Window and open `simulate_rf_pulse_b1_tool` in the `mri_rf_pulse_sim_tool`  folder, then press `F5`


## Usage:
Start by either double-clicking _simulate_rf_pulse_b1_tool.mlapp_ (should start appdesigner) or by typing appdesigner in Matlab's command window and opening the _simulate_rf_pulse_b1_tool.mlapp. Run the tool by hitting function key F5 or pressing the green "Play" button.
### Simulate RF Pulse:
3. Load an RF Pulse (either .mat or .dat in the `RF Pulses` subfolder)
    1. slr_sharp1_1ms_2000Hz will give you a Shinnar-LeRouge Pulse of 1ms and 2000Hz resulting bandwidth with "sharpness 1".
    2. AFP_2_3182ms_Hsn is an Adiabatic Full Passage Pulse of 2.3162ms (taken from Robin de Graaf's pulsewizard toolbox)
4. Press load to load the pulse.
5. You can now edit pulse duration, flipangle, pulse amplitude, and gradient strength. Bandwith factor BWfac and integration factor Sint are intrinsic properties of the pulse and can not be changed.
You can also change the range and the resolution (frequency **freq[Hz]** and position **FOV[cm]**) over which the simulation will be performed.
6. You can also change the `Sample` properties, such as **T1[s]** and **T2[s]**. You can also simulate hyperpolarized experiments by setting the polarization HP to higher than 1.
6. Pressing `Simulate!` will run the Bloch Simulation with the selected parameters.
- Note: If the parameter **yGradient Strength** is >0, a "slice selective Gradient" and a rephasing gradient (50% area of slice selective gradient) will be simulated.
7. The simulation results will be plotted in the lower left windows.
8. You can select one or multiple (hold Shift key) properties (Mx, My, ...) of the simulated magnetization.
- If you didn't load a B1 map, just use the homogeneous (homo. B1) options.
9. You can select multiple frequencies by choosing the `Show 2 Freqs`. Then you can select Freq 1 and Freq 2.
10. You can also use the slider to select different frequencies.
11. On the Bottom you can find a filed saying "save". For now I'd recommend to select all options on the list and hit the `save` button.
### Simulate repeated excitation with RF Pulse:
You can simulate repeated excitations in the Tab `Simulate Repeated Excitation (FISP/bSSFP)`. Recommended steps:
In the Tab `1D`:
1. Since you are mainly interested in the frequency response profile, set *points* of the **FOV[cm]** parameter to 1 and the *value* to 0.
2. Set the simulated frequencies (freq [Hz]) to a proper resultion (~1-5Hz).
In the Tab `Simulate Repeated Excitation (FISP/bSSFP)`:
1. Set **TR[ms]** to your desired repetition time (Note: can not be shorter than pulse duration).
2. Set number of time points of each pulse that are simulated (**N[points]**) to a reasonably high number (~500).
3. Set the number of repetitions that you want to simulate (**N[Reps]**).
4. Choose the phase difference of the RF pulse train (inc. **phase[°]**). This is typically 180 (alternating phase) or 0 (same phase).
5. Set the amplitude of the first RF pulse (First pulse amp). This should typically be 0.5.
6. Set the frequency of the RF Pulse (**f[Hz]**). This is interesting in case you want to simulate spectrally selective RF pulses. Note that for high frequency offset, you will need to increase the number of time points simulated (N[points]).
7. Set the first TR (1st **TR[ms]**). This is the duration between the first (often $\alpha$/2 RF pulse) and the second RF pulse.
- In case you choose 180° phase increment (inc. **phase[°]**), it is recommended to set the first TR to 0.5 * TR to reduced the transient phase oscillations.
- In case you choose 0° phase increment (inc. **phase[°]**), it is recommended to set the first TR to 1.0 * TR to reduced the transient phase oscillations.
8. You can choose if you want to simulate (after **N[Reps]**) 
- a "Spoiler Gradient" that sets the $M_{xy}$ component to 0.
- And/or a "Tipback pulse" tries to tipback the magnetization with a $\alpha$/2 pulse.
9. Then hit `Simulate Reps!` at the right top of the `Simulate Repeated Excitation (FISP/bSSFP)` Tab.

In case the **N[points]** parameter turns red, type in a slightly bigger or smaller number. In general, even numbers work better.This will be fixed in a future release.

You can have a look at what the frequency response profile looks like in the lower right figure. You can change the repetition with the *rep* slider.

<img src="./figures/Screenshot bSSFP Simulation example full GUI.png" alt = "Screnshot of the GUI with exemplarily parameter selection." width="1000"/>




