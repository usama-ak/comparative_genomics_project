#!/usr/bin/bash

awk '/^>/ { 
    if (length(protein) > 0) { 
        print length(protein); 
    }
    protein = ""; 
    next; 
}
{ 
    protein = protein $0; 
}
END { 
    if (length(protein) > 0) { 
        print length(protein); 
    }
}' Medicago_truncatula.MedtrA17_4.0.pep.all.fa > protein_length;
sed '1d' Medicago_truncatula_liste > new_list
paste new_list protein_length > Medicago_truncatula_protein_length;
rm protein_length;
rm new_list;
