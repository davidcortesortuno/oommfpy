import setuptools

with open('README.md') as f:
    long_description = f.read()

setuptools.setup(
    name='oommfpy',
    version='0.1',
    description=('Python lib to read and process OOMMF data files'),
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='D. Cortes',
    author_email='d.cortes@soton.ac.uk',
    packages=setuptools.find_packages(),
    install_requires=['matplotlib',
                      'numpy',
                      'click',
                      'pyvtk'
                      ],
    classifiers=['License :: BSD2 License',
                 'Programming Language :: Python :: 3 :: Only',
                 ],
    entry_points='''
        [console_scripts]
        plot_omf=oommfpy.tools.plot_slices:plot_omf_slices_cli
        omf2vtk=oommfpy.tools.omf2vtk:omf2vtk_cli
    ''',
)
