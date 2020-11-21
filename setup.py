from setuptools import find_packages, setup

setup(
	name = "networkExpansionPy",
	version = '0.0',
	description = 'Python package containing algorithms for simulating the expansion of metabolic networks',
	url = 'https://github.com/jgoldford/networkExpansionPy/',
	author = 'Joshua E. Goldford',
	author_email = 'goldford.joshua@gmail.com',
	packages = find_packages(),
	install_requires = [
		'scipy==1.1.0',
		'numpy==1.19.4',
		'pandas==0.23.4',
		'ray==1.0.1.post1'],
	include_package_data = True,
)