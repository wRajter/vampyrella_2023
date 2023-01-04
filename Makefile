# ----------------------------------
#          INSTALL & TEST
# ----------------------------------
install_requirements:
	@pip install -r requirements.txt

install_raxml-ng:
	@git clone https://github.com/amkozlov/raxml-ng.git
	@sudo apt-get install flex bison libgmp3-dev
	@cd raxml-ng
	@mkdir build && cd build
	@cmake ..
	@make

install_epa-ng:
	@sudo apt-get install autotools-dev libtool flex bison cmake automake autoconf
	@git colone https://github.com/Pbdas/epa-ng.git
	@make

install_gappa:
	@git clone --recursive https://github.com/lczech/gappa.git
	@cd gappa
	@make

check_code:
	@flake8 scripts/* ciliate_eukbank2022/*.py

black:
	@black scripts/* ciliate_eukbank2022/*.py

test:
	@coverage run -m pytest tests/*.py
	@coverage report -m --omit="${VIRTUAL_ENV}/lib/python*"

ftest:
	@Write me

clean:
	@rm -f */version.txt
	@rm -f .coverage
	@rm -fr */__pycache__ */*.pyc __pycache__
	@rm -fr build dist
	@rm -fr ciliate_eukbank2022-*.dist-info
	@rm -fr ciliate_eukbank2022.egg-info

install:
	@pip install . -U

all: clean install test black check_code


# ----------------------------------
#      UPLOAD PACKAGE TO PYPI
# ----------------------------------
PYPI_USERNAME=<AUTHOR>
build:
	@python setup.py sdist bdist_wheel

pypi_test:
	@twine upload -r testpypi dist/* -u $(PYPI_USERNAME)

pypi:
	@twine upload dist/* -u $(PYPI_USERNAME)