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
    
        - :fontawesome-solid-exchange-alt: - Recalculated high confidence miRNA set: 2162 (version 22) :octicons-arrow-right-16: 3320 (version 22.1)
        
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
        
        [Detailed change list :material-update:](v22-diff.md){: target="_blank" .md-button .md-button--primary }           
            
    === "21"
       
        **Overview**:
            
        - Version released in June 2014
        - 28645 entries representing miRNA precursor 
        - 35828 entries representing mature miRNA 
        - 223 species in which all entries are found 
                
        **Differences from previous version**:
        
        - :material-plus: - 4196 new precursor sequences
        - :material-plus: - 5441 new mature products
        - :material-plus:  - 17 new species
        - :material-plus: - High confidence miRNA set was created (1996 entries labeled as of high confidence)    
        - :material-rename-box: - 169 precursors and 353 mature miRNAs have changed names
        - :material-delete-forever: - 72 misannotated and duplicate sequences have been deleted                     
        
        !!! warning "Issue"
            Some of .gff files are different comparing to version 22 or above. There are possible inconsistencies and differences in chromosome names.  
        
        [Detailed change list :material-update:](v21-diff.md){: target="_blank" .md-button .md-button--primary }
        
    === "20"
        
        **Overview**:
            
        - Version released in June 2013
        - 24521 entries representing miRNA precursor 
        - 30424 entries representing mature miRNA 
        - 206 species in which all entries are found 
                
        **Differences from previous version**:
        
        - :material-plus: - 3355 new precursor sequences
        - :material-plus: - 5393 new mature products
        - :material-plus:  - 13 new species   
        - :material-rename-box: - 144 precursors and 810 mature miRNAs have changed names
        - :material-delete-forever: - 98 misannotated and duplicate sequences have been deleted                     
        
        !!! warning "Issue"
            Some of .gff files are different comparing to version 22 or above. There are possible inconsistencies and differences in chromosome names.
        
        !!! missing "Missing feature"
            This version does not define high confidence set of miRNA, thus `high_confidence` fields will be set to default value which is `False`.               
        
        [Detailed change list :material-update:](v20-diff.md){: target="_blank" .md-button .md-button--primary }
        
    === "19"
           
        **Overview**:
            
        - Version released in August 2012
        - 21264 entries representing miRNA precursor 
        - 25141 entries representing mature miRNA 
        - 193 species in which all entries are found 
                
        **Differences from previous version**:
        
        - :material-plus: - 3171 new precursor sequences
        - :material-plus: - 3625 new mature products
        - :material-plus:  - 25 new species   
        - :material-rename-box: - 313 precursors and 4666 mature miRNAs have changed names
        - :material-rename-box: - Mature sequences from all species are now designated -5p and -3p, rather than miR/miR*.
        - :material-delete-forever: - over 130 misannotated and duplicate sequences have been deleted                     
        
        !!! warning "Issue"
            Some of .gff files are different comparing to version 22 or above. There are possible inconsistencies and differences in chromosome names.
    
        !!! warning "Issue"
            Since version 19, the files with organisms data are lacking NCBI taxonomy ID. Without this information function `get_taxid` would be rendered obsolete and data structure redesign should be conducted. To keep the intended functionality and structural consistency alongside further compatibility, current taxonomy IDs are retrieved from NCBI database based on organism names availible in data files from version 19. This largely solves the problem, but some organism names were changed and as a result some records might be lacking in organism/taxonomy data.           
        
        !!! missing "Missing feature"
            This version does not define high confidence set of miRNA, thus `high_confidence` fields will be set to default value which is `False`.               
        
        [Detailed change list :material-update:](v19-diff.md){: target="_blank" .md-button .md-button--primary }
        
    === "18"
            
        **Overview**:
            
        - Version released in November 2011
        - 18226 entries representing miRNA precursor 
        - 21643 entries representing mature miRNA 
        - 168 species in which all entries are found 
                
        **Differences from previous version**:
        
        - :material-plus: - 1488 new precursor sequences
        - :material-plus: - 1929 new mature products
        - :material-rename-box: - 50 precursors and 1417 mature miRNAs have changed names
        - :material-rename-box: - Mature sequences from all human, mouse, and C. elegans precursors are now designated -5p and -3p, rather than miR/miR*.
        - :material-delete-forever: - over 45 misannotated and duplicate sequences have been deleted                     
        
        !!! warning "Issue"
            Some of .gff files are different comparing to version 22 or above. There are possible inconsistencies and differences in chromosome names.
    
        !!! warning "Issue"
            Since version 19, the files with organisms data are lacking NCBI taxonomy ID. Without this information function `get_taxid` would be rendered obsolete and data structure redesign should be conducted. To keep the intended functionality and structural consistency alongside further compatibility, current taxonomy IDs are retrieved from NCBI database based on organism names availible in data files from version 19. This largely solves the problem, but some organism names were changed and as a result some records might be lacking in organism/taxonomy data.
    
        !!! warning "Issue"
            Mir-Us sometimes relies on names and different naming convention is suspected to render some records lacking names or will cause assignment of the incorrrect name from the perspective of newer versions (e.g. certain names can be incorrectly trimmed) to which Mir-Us is native.                       
        
        !!! missing "Missing feature"
            This version does not define high confidence set of miRNA, thus `high_confidence` fields will be set to default value which is `False`.               
        
        [Detailed change list :material-update:](v18-diff.md){: target="_blank" .md-button .md-button--primary }