from setuptools import setup, Extension

setup(
    name = 'makegms',
    version = '1.5.5',
    description = 'Fast and easy tool for genome mappability estimation',
    url = "https://github.com/evgeny-bakin/GeMaTrIA",
    author = 'Igor V.',
    author_email = 'latur@me.com',
    ext_modules = [Extension(
      'makegms', 
      ['src/makegms.c'], 
      include_dirs = ['src'],
      extra_compile_args=['-std=c99', '-m64', '-O', '-lpthread']
    )],
    py_modules = ['makegms']
)

# rm -rf makegms* dist build
# python3 setup.py build
# python3 setup.py install
