import setuptools

with open("README.md",'r',encoding='utf-8') as fh:
    long_description=fh.read()


setuptools.setup(
    author ="Zhetao Tan; Xinyi Song; Lijing Cheng; Rebecca Cowley, Huifeng Yuan, Guilherme Castelao, Simona Simoncelli, Shoichi Kizu, Ricardo Locarnini, Tim Boyer, Franco Reseghetti, Viktor Gouretski",
    long_description = long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    package_data={'DC_OCEAN':['util/*','support/*.py','tests/Examples_netCDF_files/*','Input_files/WOD18_sample_1995/*','Input_files/*.txt','tests/*.py']},
    packages=setuptools.find_packages(),
    install_requires=['numpy >= 1.19.1','netCDF4 >= 1.5.4','timezonefinder >=6.0.1','pandas >=1.0.3','scipy >=1.7.3','argparse >=1.4.0'],
)
