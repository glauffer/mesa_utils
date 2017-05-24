from setuptools import setup

setup(name='mesa_utils',
      version='0.2',
      description='Functions to work and plot MESA data',
      url='http://github.com/glauffer/mesa_utils',
      author='Gabriel Lauffer Ramos',
      author_email='gabriellramos@gmail.com',
      license='MIT',
      packages=['mesa_utils'],
      install_requires=[
          'numpy',
          'matplotlib'
          'mesa_reader'
      ],
      dependency_links=['https://github.com/wmwolf/py_mesa_reader'],
      zip_safe=False,)
