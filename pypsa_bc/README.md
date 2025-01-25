# BC_Combined_Modelling
Project combining end-to-end workflow combining BC_Nexus, PyPSA_BC and the Bi-directional Linking Tool.

The standalone components (gits) of the projects are -
- [BC_Nexus](https://github.com/DeltaE/BC-CLEWS-Model).
- [BC_PyPSA](https://github.com/DeltaE/PyPSA_BC).
- [Linking Tool](https://github.com/DeltaE/Linking_tool).


## Developers

- **Bruno Borba**  
  - Concept build and code development for storage and hydro modelling, timeseries aggregation techniques.
  - input data preparation automation for CLEWs model.

- **Elias Islam**  
  - Combined Modelling workflow and coding structure
  - Model linking tool
  - Resource capacity disaggregation tool  
  
- **Pierre McWhannel**  
  - PyPSA model and workflow development
  - Storage and Hydro modelling conceptualization for CLEWs model.
  - Timeslice aggregation conceptualization for CLEWs model.

## Acknowledgements

- **Dr. Taco Niet**  
  - Principal Investigator (PI) of the project  
  - Overall supervision of the development  
  - Project development and funding arrangements from PICS

## Affiliation

- Delta E+ Research Lab, Simon Fraser University (SFU)

## Version

- **No:** 1  
- **Release:** 202409

Developer_Remarks:
 Key_upgrade:
    - Input data preparation automation add-on to CLEWs framework (tailored for BCNexus Model).
 New_components:
   - Battery and Hydro Storage, Timeslice up/down scaling upto 1-hour resolution, cascaded hydro power.
"""
