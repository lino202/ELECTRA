# ELECTRA Changelog
===================
Record of changes for ELECTRA project.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
ELECTRA uses a vMAJOR.MINOR.PATCH versioning scheme. The
* MAJOR version is increased when functionality that is not backwards compatible is added, the
* MINOR version is increased when backwards compatible functionality is added, and the
* PATCH version is increased when bug fixes are added.

## Version History
------------------
* v0.6.3
* v0.6.2
* v0.6.1
* v0.6.0
* v0.5.2
* v0.5.1
* v0.5.0
* v0.4.5
* v0.4.4
* v0.4.3
* v0.4.2
* v0.4.1
* v0.4.0
* v0.3.1
* v0.3.0
* v0.2.0
* v0.1.0

## 2025-01-02: v0.6.3
---------------------

### Added
* Now checkpointing (load and save states) is available, still this might not function when load_curve is used for dynamically modifying a parameter
* We added the paci2020 cellular model (10.1016/j.bpj.2020.03.018), please pay attention to the model comments in paci2020.cpp and new comments in paci_ventri.cpp (corresponds to paci model 2013), specially regarding the stim amplitude and prepacing.
* Now ElectraCell recieves the dt parameter which is the integration step. The output dt in ElectraCell is always 0.01 ms.
* New examples for depicting new features and the use of purkinje.


### Fixed
* States and params as well as valuesin the cell models where checked for the expression with decimals (1 is 1.) so the calculations did not yield rounded int values (1/1000 = 0 but 1./1000. = 0.001)
* Stewart model used the exp function from eigen with yielded equal results as std::exp, but now std::exp is used for consistency with other models.


### Notes
* For consistency new cell models should be added at the end of the enum list EpModelType (ep_basic.hpp) as checkpointing use this index, if a version of electra with a new cellular model (added in the middle of this list) is used for LOADING the states of a cellular model which was at the end of the list -> the index has changed and the wrong cellular model will be loaded which will yield errors and the necessity to change the model index in the .elc files for compability with the new electra

## 2024-07-13: v0.6.2
---------------------

### Added
* The Ohara model is fixed. Now it stabilizes in Electra as well, we fixed the sign of Istim in dki. Also we fixed the gKb for epi, the tau_rels and the nca integration method. Now Electra has the default Ohara cell model "ohara2011" and "ohara2011_inatt" which defines the INa as done by Ten-tusscher (recommended by Ohara in the notes of his paper).
* The Stewart model was double checked and left as the default as state variables initial values were different and ordering of the calculations. We fixed the definitions of prm and var which saves memory specially in fine meshes.
* The GNa in the Gaur2021 model has been left as default as the depolarization and propagation seem normal in tissue. We left the model as before,
with only the ventricular default cell type, if you want to make changes you can use manual init file in ElectraCell and ElectraSim over the variables that are already define, you can obsiously change the model and compile the code but for now the default variables are these.
* Other minimal fixes and additions, see commits.

## 2024-03-14: v0.6.1
---------------------

### Added
* SingleCellDev app now is called ElectraCell, and now is available again
* ElectraCell now is an app that receives arguments as input for simulating single cells, the arguments format is 
./ElectraCell /path/file_name ep_model cell_type sim_time stim_start stim_dur stim_cycle_length stim_amp
* Now ElectraCell saves the manual init file for using with ElectraSim
* The GNa in the Gaur2021 model has been double as Electra simulates mainly tissue

## 2024-02-26: v0.6.0
---------------------
### Notes
* In general, ElectraSim now is more suitable in terms of memory for working with bigger meshes and longer periods of time. ElectraSim now saves states dynamically in binary ensight gold format. 
* The code and code structure now is a little bit simpler and it has a more clear way of building and using it in release and debug modes.

### Added
* Read fibers as a separated .json file
* Save ensight binary (saves dynamic (RAM) and storage (ROM) memory)
* Dynamic saving during simulation computation of states (save dynamic memory)
* BLOCK_CELL_CURRS for generating currents blocks or not. Now, current blocks are desativated as default (saves dynamic memory). Also majority of currents have a conductivity associated to them. Pay attention that in the tentusscher2006 and Stewart2009 models the current block can be set but it does not actually perform the block as it is not defined in the source code.

### Fixed
* Minor fixes in computation of CS diffusivities in the end branches
* Setting a more clear instruction for building ElectraSim
* TinyXML is not need in Electra anymore as the paraview exporter is deprecated but it is still in CLOUDEA. But we do not install tinyxml anymore we basically added its .c and .h files and complied CLOUDEA with it.

### Deprecated/Desactivated
* In general, several situations are now deprecated/desactivated majorly because we change the way we save the states for saving memory. 
* Desactivate: The bidomain model as the Compute member function needs to implement the dynamic saving as for the monodomain setting. 
* Desactivate: The MCM as the FictitiousValuesToReal function requires the transmembrane potential for all time steps and nodes. As we noe save dynamically this values has to be read again after simulating which is not implemented for now but it is not difficult to do.
* Desactivate: ElectraPre was desactivated for simplifyng the building but it can be uncommented in the CMakeList.txt file under /apps. ElectraPost was desactivated as well as it is just as skeleton for now.
* Deprecated: The postprocessing, as APD was not working, we need to read back all states and we wanted to isolate ElectraSim for only simulating. Postprocessing can be done (ans is majorly done) with thirdparty/in-house made code or by the app ElectraPro in the future.
* Deprecated: the paraview and ascii exporters. There is no need to save the data in ascii. All is binary now.


## 2022-02-20: v0.5.2
---------------------
### Fixed
* Corrected diffusion transition nodes of conduction system. Transition branches were 10 nodes long before. Now they are 10 nodes long except if a bifurcation is
  found. In this case the last transition node is the one located at the bifurcation.



## 2022-02-07: v0.5.1
---------------------
### Fixed
* Corrected an issue with the determination of fiber direction unit vectors. Now Electra normalizes the given fiber direction vectors internally


## 2022-01-02: v0.5.0
---------------------
### Added
* New applications: ElectraSim and ElectraPre

### Added
* Old application: ELECTRA-console

### Changed
* Json script interface. Check tutorial and examples folders for more info


## 2021-03-04: v0.4.5
---------------------
### Added
* Maleckar2009 electrophysiology model (Mathematical simulations of ligand-gated and cell-type specific effects on the action potential of human atrium, Maleckar et al. 2008, Progress in biophysics and molecular biology)

### Changed
* Included modifications on the IKCa current of the Courtemanche1998 model according to Engel and Peñaranda provided by Chiara Celotto

### Fixed
* Improved integration of gate variable equations in Courtemanche1998 model
* Improved integration of gate variable equations in Bueno2008 model
* Improved integration of gate variable equations in Stewart2009 model
* Adaptive diffusion timestep could exceed the maximum timestep in some cases. Now this is corrected.


## 2021-01-21: v0.4.4
---------------------
### Added
* Modified version of the Gong2020 electrophysiology model using the gate variable equations for INa from Tentusscher2006 model

### Fixed
* Improved integration of gate variable equations in Grandi2011a model


## 2021-01-20: v0.4.3
---------------------
### Added
* Tentusscher2006 electrophysiology model (Alternans and spiral breakup in a human ventricular tissue model , Tentusscher et al. 2006,  Am J Physiol Heart Circ Physiol)

### Fixed
* Improved integration of gate variable equations in Gong2020 model using the Rush Larsen method


## 2021-01-01: v0.4.2
---------------------
### Added
* Gong2020 electrophysiology model (Quantitative analysis of variability in an integrated model of human ventricular electrophysiology and b-adrenergic signaling, Gong et al. 2020, JMCC)

### Changed
* Objects *ap models* and *ap type* in input script (.json) are renamed to *ep models* and *ep type*
* electrophysiology model *Ohara* is renamed to *Ohara2011m*
* electrophysiology model *Bueno* is renamed to *Bueno2008*
* electrophysiology model *PaciVentri* is renamed to *Paci2013v*
* electrophysiology model *Courtemanche* is renamed to *Courtemanche1998*
* electrophysiology model *GrandiAtri* is renamed to *Grandi2011a*
* electrophysiology model *Maccannell* is renamed to *Maccannell2007*
* electrophysiology model *Stewart* is renamed to *Stewart2009*


## 2020-07-19: v0.4.1
---------------------
### Added
* Assignment of conduction system tree diffusivity on nodesets


## 2020-06-16: v0.4.0
---------------------
### Added
* Binary output of cells state at the end of the simulation
* Ability to start the simulation with the cells state initialized by loading a binary file
* Implementation of the ventricular conduction system

### Changed
* Significant changes on the Electra interface. Now there are two sections for tissue and conduction system setup
* Changed the conductivity attribute in the material section to diffusivity

### Fixed
* Corrections in meshfree implementation
* Corrected fibers assignment for homogeneous fiber orientation


## 2020-06-08: v0.3.1
---------------------
### Added
* Paci et al 2013 action potential model for stem cell derived human ventricular cell

### Changed
* The sign of the stimulus amplitude is now reversed. A positive amplitude value is required to trigger cell depolarization

### Fixed
* Corrections in meshfree implementation
* Corrected fibers assignment for homogeneous fiber orientation


## 2020-05-19: v0.3.0
---------------------
### Added
* Support for immersed grid Mixed Collocataion meshfree models
* Support for different dilatation coefficient for surface and interior nodes for meshfree support domains
* Calculation of APD in post processing

### Fixed
* Correction in cAF parameter application in Courtemanche action potential model


## 2020-03-19: v0.2.0
---------------------
### Added
* Time-varying action potential model parameters using a loading curve

### Changed
* Replacement of [Armadillo](http://arma.sourceforge.net/) library dependency with [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)


## 2020-02-16: v0.1.0
---------------------
* Experimental private release
