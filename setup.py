import setuptools

with open("README.md",'r',encoding='utf-8') as fh:
    long_description=fh.read()


setuptools.setup(
    name="DC_OCEAN",
    version="1.3.4",
    author ="Zhetao Tan; Xinyi Song, Lijing Cheng, Rebecca Cowley, Huifeng Yuan, Guilherme Castelao, Simona Simoncelli, Shoichi Kizu, Ricardo Locarnini, Tim Boyer, Franco Reseghetti, Viktor Gouretski",
    author_email = "tanzhetao19@mails.ucas.ac.cn; songxinyi231@mails.ucas.ac.cn",
    description = "DC_OCEAN: An algorithm to detect the ocean in-situ duplicate profiles (Song et al., 2024, FMS)",
    long_description = long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/IQuOD/duplicated_checking_IQuOD/",
    include_package_data=True,
    package_data={'DC_OCEAN':['util/*','support/*.py','tests/Examples_netCDF_files/*','Input_files/WOD18_sample_1995/*','Input_files/*.txt','tests/*.py']},
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    install_requires=['numpy >= 1.19.1','netCDF4 >= 1.5.4','timezonefinder >=6.0.1','pandas >=1.0.3','scipy >=1.7.3','argparse >=1.4.0'],
    python_requires='>=3.8',
)
