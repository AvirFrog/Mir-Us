# Mir-Us
> This is a placeholder of a README. Full version will be released with the 1.0 version of the Mir-Us.

 Fast, intuitive mirBase library, for solving problems with mirBase

## Getting started
To get Mir-Us to work, it is necessary to import `miBase` package and then create a database object:
```python
import miBase
m = miBase.MiRBase()
```
miRBase version can be specified while creating a new object:
```python
import miBase
m = miBase.MiRBase(version="20")
```
Example of a full usage:
```python
sample_gallus = m.get_mirna(mirna_id=["MIMAT0001185", "MIMAT0025825", "MIMAT0007451"])
print(sample_gallus)
```
More information about usage and technical aspects of the Mir-Us can be acquired in the docs page.
Documentation includes:
- getting started guide
- user cookbook
- full compatibility description between versions of miRBase
- technical reference documentation from code

Docs page is not hosted externally yet, thus it must be initialised locally ([instructions](#accessing-docs)).

## Accessing docs
### Installing dependencies
To access docs, a MkDocs, material theme for MkDocs and mkdocstrings package must be installed.

#### MkDocs installation:
```commandline
$ pip install mkdocs
```
#### material theme for MkDocs installation:
```commandline
$ pip install mkdocs-material
```
#### mkdocstrings installation:
```commandline
$ pip install mkdocstrings
```
### Running docs
To run docs, you have to be in the main directory of the Mir-Us and initialise the MkDocs. Then head to the localhost (http://127.0.0.1:8000/) and documentation should appear in your default web browser.
```commandline
$ mkdocs serve
```

## Authors

- **Dudczak Kacper**
- **Michalczyk Maciej**

## Credits
- **Dudczak Kacper**
- **Michalczyk Maciej**
- **Marta Wysocka**
- **Marek Żywicki**

