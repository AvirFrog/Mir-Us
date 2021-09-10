# Getting started
This section shows an example usage of the Mir-Us package for the absolute beginner. There are a few steps that have to be done to receive data from [miRBase](https://www.mirbase.org/).

## Import package
To use the package, it must be first imported:
> to się może zmienić
    ```python
    from mirus import miBase
    ```
## Initialize miBase object
In order to call functions, a miBase object must be initialized:

!!! example "Example of basic initialization"
    ```python
    m = miBase.MiRBase()
    ```

This code above will create a miBase object which will provide functions to retrieve data from the current version of miRBase. However, one might be working on an old project, which requires access to older versions of miRBase. It is possible with Mir-Us.

!!! info
    Mir-Us supports older versions of miRBase. To use them, the `version` parameter must be provided. Backwards compatibility does not end on version 20, however. [Details :octicons-link-16:](versions.md){ .md-button .md-button--primary }

!!! example "Example of initialization with version defining"
    ```python
    m = miBase.MiRBase(version="20")  # for object 'm' the data from version 20 will be used.
    ```

!!! warning "The lower the version number, the older database is and it will provide fewer data. In extreme cases, incomplete data might be returned."

## Receive data
Receiving data requires the usage of provided functions, each of them is designed for a specific task.

For example, one might want to obtain information about existing miRNAs on a human's chromosome X upstream from 153300000th nucleotide. For this task, the function `get_mirna` is the most suitable.

!!! example "Example of miRNA search on human's chromosome X from specific genomic location"
    === "Code"
        ```python
        data = m.get_mirna(chr="chrX", organism_name="Homo sapiens", start="153300000")
        ```
    === "Result"
        ```
        [Mir-Us]  'get_mirna' found 9 results in 0.091857 seconds
        ```

In the 'Result' tab a search report was printed to stdout and most importantly, a list of miRNA objects was returned, which contain all the parsed information.

## Next steps
With this basic knowledge it is now possible to learn Mir-Us in-depth. If one is looking for specific application of this tool, then it is suggested to read the 'User cookbook'. If one is looking for more specific information about some functionality, then it is advisable to read the 'Reference documentation'.

**It's time to choose!**

[User cookbook :fontawesome-solid-book-medical:](cookbook.md){ .md-button .md-button--primary }   [Reference documentation :octicons-book-16:](miObject.md){ .md-button .md-button--primary }