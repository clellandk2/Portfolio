from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / "README.md").read_text(encoding="utf-8")


setup(
    name="VISION_Multiphonon_Correction",  
    version="1.0.0",  
    description="A code to perform multiphonon corrections on neutron scattering data from an indirect geometry instrument.",  
    long_description=long_description,  
    long_description_content_type="text/markdown",  
    url="https://github.com/clellandk2/Portfolio/tree/main/Multiphonon%20Correction",  
    author="Kevin Clelland",  
    keywords="Neutron_Scattering, Indirect_Geometry", 

    packages=find_packages(exclude = "Example_Data"),
    install_requires=[
        "numpy",
        "matplotlib",
    ],

    python_requires=">=3.7, <4",    
        entry_points={  
        "console_scripts": [
            "sample=VISION_Multiphonon_Correction.VISION_Multiphonon_Correction:setupSQEs",
        ],
    },
)