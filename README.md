# ThermalDex
[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIt)
![Qt](https://img.shields.io/badge/Qt-%23217346.svg?style=for-the-badge&logo=Qt&logoColor=white)
![HotCookie](./HouseKeeping/GitHubBadges/HotCookie.svg)

ThermalDex is an application developed in Python using RDKit and PyQt5 that facilitates the identification and assessment of thermal hazards in molecules. It offers suite of tools for analyzing chemical structures, calculating safety metrics, and generating detailed reports to aid in the safe handling and processing of potentially hazardous materials.

## Features

- **SMILES Input**: ThermalDex accepts input in the form of SMILES (Simplified Molecular Input Line Entry System), allowing users to input molecular structures easily. ThermalDex will also accept direct copy/paste from Chemaxon's Marvin/JChem files (ChemDraw Support is *ToDo*).

- **High Energy and Explosive Group Identification**: The application identifies high-energy and potentially explosive groups within molecules, providing valuable insights into their thermal stability and safety. This is achived via a substructure search and a .csv of the relavent groups.

- **O.R.E.O.S. Assessment**: ThermalDex calculates O.R.E.O.S. assessment values,<sup>1</sup> incorporating safety metrics such as Oxygen Balance and Rule of Six to evaluate the potential hazards associated with a molecule.

- **DSC Data Analysis**: The application can record and analyze the results from Differential Scanning Calorimetry (DSC), enabling users to determine critical parameters such as T<sub>D24</sub> (Temperature of Decomposition),<sup>2</sup> Yoshida Impact Sensitivity, and Yoshida Explosive Propagation.<sup>3</sup>

- **Graphical Visualization**: ThermalDex features a Qt GUI that displays molecular structures graphically, and the calculated hazard levels are color-coded for intuitive interpretation.

- **PDF Report Generation**: The application generates detailed PDF reports of thermal hazard assessments for each molecule, providing comprehensive documentation of the analysis results.

- **Local Storage**: ThermalDex maintains a user-friendly local .csv file that stores all assessment results, making it easy to track and manage data over time. ThermalDex uses sqlite to store all results and changes to help prevent data loss. This can be found as 'ThermalDex.db' in the _core folder.

<p><small>[1]: <i>Org. Proc. Res. Dev.,</i> 2021, 25, 2, 212-224</small><br>
<small>[2]: <i>Angew. Chem. Int. Ed.,</i> 2020, 59, 15798-15802</small><br>
<small>[3]: <i>Org. Proc. Res. Dev.,</i> 2020, 24, 1, 67-84</small></p>

## Getting Started

### Executable Releases 
Executables will be released via Github Releases (currently we are testing internally but will release in near future). These executables have be generated thanks to Pyinstaller. ThermalDex has been a Windows first project due to the nature of the devices used for testing, but I will get a MacOS executables shortly and aim to generate at least a .deb package for Linux down the line. If compatabliity issues are encountered I would recomemend looking to running in Python. 

### Python Code

To get started with ThermalDex, follow these steps:

1. Clone the repository to your local machine:
   ```
   git clone https://github.com/Tomisall/ThermalDex.git
   ```

2. Install the required dependencies using pip:
   ```
   pip install -r requirements.txt
   ```

4. Launch the application:
   ```
   python ThermalDex.py
   ```

5. Input the SMILES notation of the molecule you wish to analyze and explore the various features and tools offered by ThermalDex. Note you will need to add at least one compond before trying the search feature, otherwise this may result in a crash.

## Contributing

Contributions to ThermalDex are welcome! If you encounter any issues, have suggestions for improvements, or would like to contribute new features, please feel free to submit a pull request or open an issue on GitHub.

Please get in touch if you have any issues or would like to see any additional features. I want ThermalDex to be a really useful tool for those thinking about (and those who don't want to but should think about) Thermal Hazards in organic synthesis.

## Disclaimer

This python code was written entirely by a synthetic organic chemist who was also working in the lab ('aka his actual job') and raising his new Son (who his about 2 weeks older than this project and is adorable). It could do with being neater and having a good old refactor. My ability to name variables/methods/functions consistently is admittedly not great (Its hard to stick to snake or camel case when you use functions from different import modules which use different standards, but the blame is all mine). I will work on getting the code base tidied up so please bare with me.

## License

ThermalDex is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.