# Contact-Based binding affinity estimation using simplified surface area decomposition

Intermolecular binding affinity is a critical measure in structural biology, drug discovery, and biomolecular engineering. Among various physics-based and empirical approaches, surface area-based scoring models offer a computationally efficient and interpretable strategy. These models rely on the observation that the formation of a stable complex is associated with the burial of solvent-accessible surface area between interacting molecules.

## Theoretical Background

### Contact Surface Area

The **solvent-accessible surface area (SASA)** is the surface traced by the center of a water-sized probe (typically 1.4 Å radius) rolling over the molecular surface. The contact area between two molecular components A and B can be approximated using:

$A_{contact}(A,B) = \frac{1}{2} \left( SASA(A) + SASA(B) - SASA(A \cup B) \right)$

This equation estimates the buried surface area resulting from the interaction.

> *Note: The ½ factor accounts for symmetry in the interface, avoiding double-counting the buried area from each molecule. Some sources omit this factor when estimating total buried SASA instead of interfacial contact area.*

### Surface Type Decomposition

In this simplified model, atoms are grouped into two classes:

* **Polar (P):** Atoms likely to participate in hydrogen bonding or carry formal/partial charges (e.g., N, O, S,  P, charged side chains)
* **Nonpolar (NP):** Aliphatic and aromatic carbon atoms that are generally hydrophobic

Three interaction types are considered:

* P\:P (polar-polar)
* NP\:NP (nonpolar-nonpolar)
* P\:NP (mixed polar-nonpolar, includes both P from A and NP from B, and vice versa)

Each of these contact types is associated with an empirical weight reflecting their typical contribution to binding free energy.

### Scoring Function

The surface binding affinity estimate (in kcal/mol) is computed as:

$$
Affinity_{contact} = \alpha_{P:P} \cdot A_{P:P} + \alpha_{NP:NP} \cdot A_{NP:NP} + \alpha_{P:NP} \cdot A_{P:NP}
$$

Where:

* $A_{i:j}$: Contact surface area between atom types $i$ and $j$
* $\alpha_{i:j}$: Empirical coefficients (kcal/mol/Å²) for each type of contact

#### Default Coefficients

Based on literature and empirical scoring models:

| Interaction Type | Coefficient $\alpha_{i:j}$ (kcal/mol/Å²) | Rationale/Source                                                                                                                                        |
| ---------------- | ---------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| P\:P             | -0.015                                   | Favorable H-bonding or polar interactions ([Eisenberg & McLachlan, 1986](https://doi.org/10.1038/319199a0))                                             |
| NP\:NP           | -0.025                                   | Hydrophobic effect, burial of apolar surfaces (\[Eisenberg & McLachlan, 1986]; [Zhang et al., 1997](https://doi.org/10.1016/S0022-2836%2897%2900108-2)) |
| P\:NP            | -0.005                                   | Mixed interactions are less favorable and sometimes destabilizing (\[Zhang et al., 1997])                                                               |

These values can be modified to suit a particular force field or calibrated against experimental data.

## Implementation in VMD

The provided `contact_surface.tcl` script performs the following steps:

1. Defines selections for molecules A and B.
2. Classifies atoms in each selection into polar (P) and nonpolar (NP). You can set your own classification rules.
3. Computes surface contact areas between each interaction type per trajectory frame.
4. Applies the scoring function using the defined $\alpha$ coefficients.
5. Writes the output to a tab-delimited file with per-frame data.


## Use Cases

* Screening protein-protein, protein-membrane protein-ligand, protein-DNA... complexes for interface quality
* Monitoring binding/unbinding during MD trajectories
* Estimating relative binding affinity for mutant or modified structures

## Limitations and Considerations

* This model is qualitative and does not include entropic contributions.
* Water-mediated contacts are not explicitly considered.
* Atom classification is rule-based and may require adjustment for specific force fields.
* Coefficients $\alpha$ are general approximations; context-specific calibration is recommended for quantitative predictions.

## References

* Eisenberg, D., & McLachlan, A. D. (1986). Solvation energy in protein folding and binding. *Nature*, 319(6050), 199–203. [https://doi.org/10.1038/319199a0](https://doi.org/10.1038/319199a0)
* Zhang, C., Liu, S., Zhu, Q., & Zhou, Y. (1997). A knowledge-based energy function for protein-ligand, protein-protein, and protein-DNA complexes. *J. Mol. Biol.*, 267(3), 707–726. [https://doi.org/10.1016/S0022-2836(97)00108-2](https://doi.org/10.1016/S0022-2836%2897%2900108-2)
* Connolly, M. L. (1983). Solvent-accessible surfaces of proteins and nucleic acids. *Science*, 221(4612), 709–713.

---

