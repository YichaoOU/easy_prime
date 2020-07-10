from setuptools import setup, find_packages
import easy_prime
try:
	with open("README.md", "r") as fh:
		long_description = fh.read()
except:
	long_description = ""
setup(
	name='easy_prime',
	version=str(easy_prime.__version__),
	description="Prime editor gRNA design tool",
	author="Yichao Li",
	author_email='Yichao.Li@stjude.org',
	url='https://github.com/YichaoOU/easy_prime',
	packages=find_packages(),
	license='LICENSE',
	scripts=['bin/easy_prime','bin/easy_prime_vis'],
	package_data={'': ["*","test/*","model/*"]},
	include_package_data=True,
	long_description=long_description,
	long_description_content_type='text/markdown'	,
	keywords='prime editor',
	classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Topic :: Scientific/Engineering :: Visualization',
		'Topic :: Scientific/Engineering :: Information Analysis',
		'Operating System :: Unix',
		'Natural Language :: English',
		"Programming Language :: Python :: 3"
	]
)

