# import os
from distutils.command.build_ext import build_ext
from distutils.core import Extension

# See if Cython is installed
try:
    from Cython.Build import cythonize
# Do nothing if Cython is not available
except ImportError:
    # Got to provide this function. Otherwise, poetry will fail
    def build(setup_kwargs):
        print('Not compiling Cython module')
        pass
# Cython is installed. Compile
else:
    # This function will be executed in setup.py:
    def build(setup_kwargs):
        # The file you want to compile
        com_args = ['-std=c99', '-O3']
        extensions = [
            Extension("oommfpy.tools.clib",
                      ["oommfpy/tools/clib/clib.pyx",
                       "oommfpy/tools/clib/vtk_writer.c"],
                      extra_compile_args=com_args
                      ),
        ]

        # gcc arguments hack: enable optimizations
        # os.environ['CFLAGS'] = '-O3'

        # Build
        setup_kwargs.update({
            'ext_modules': cythonize(
                extensions,
                language_level=3,
                compiler_directives={'linetrace': True},
            ),
            'cmdclass': {'build_ext': build_ext}
        })
