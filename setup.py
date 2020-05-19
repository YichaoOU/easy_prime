from setuptools import setup, find_packages

with open("README.md", "r") as fh:
	long_description = fh.read()

setup(
	name='easy_prime',
	version='1.0',
	description="Prime editor gRNA design tool",
	author="Yichao Li",
	author_email='Yichao.Li@stjude.org',
	url='https://github.com/YichaoOU/easy_prime',
	packages=find_packages(),
	license='LICENSE',
	scripts=['bin/easy_prime'],
	package_data={'': ["*","test/*","model/*"]},
	include_package_data=True,
	long_description=long_description,
	long_description_content_type='text/markdown'	,
	keywords='prime editor',

)

