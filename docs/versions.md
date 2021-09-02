# Compatibility
While designing Mir-Us, we knew that it can't be constrained to the newest version only. We are aware, that many researchers are still using data from a few years ago. That way, we are reaching these people with backwards compatibility that spans as long as 10 years back! But it must be kept in mind, that such vast compatibility is not possible without some trade-offs.

!!! info "Supported versions and their features"

    === "22.1 ('CURRENT')"
    
        **Overview**:
        
        - The newest version (December 2018)
        - The amount of entries is identical with version 22:    
            - 38589 entries representing miRNA precursor 
            - 48885 entries representing mature miRNA 
            - 271 species in which all entries are found
            
        **Differences from previous version**:
    
        - Recalculated high confidence miRNA set: 2162 (version 22) :octicons-arrow-right-16: 3320 (version 22.1)
        
        !!! check "No known issues"            
            
    === "22"
    
        **Overview**:
            
        - Version released in March 2018
        - 38589 entries representing miRNA precursor 
        - 48885 entries representing mature miRNA 
        - 271 species in which all entries are found 
                
        **Differences from previous version**:
        
        - :material-plus: - 10031 new precursor sequences
        - :material-plus: - 13149 new mature products
        - :material-plus:  - 48 new species
        - :material-rename-box: - 115 precursors and 496 mature miRNAs have changed names
        - :material-delete-forever: - 87 misannotated and duplicate sequences have been deleted                     
        - :fontawesome-solid-exchange-alt: - Recalculated high confidence miRNA set: 1996 (version 21) :octicons-arrow-right-16: 2162 (version 22)
        
        !!! check "No known issues"   
        
        [Detailed change list :material-update:](v22-diff.md){ .md-button .md-button--primary }           
            
    === "21"
    === "20"
    === "19"
    === "18"