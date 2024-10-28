# FAIR Chemistry Data

## Research Question

Can we predict the activity of a compound against a kinase?

## Data Types

### Essential

* Compound InChI (string) <-- this is our unique identifier
* Compound Smiles (string)
* Compound InChI (string)
* List of compound synonyms (and their type) (list of lists containing strings, e.g. [[name 1, research code], [name 2, IUPAC Name]])
* Protein Target (list of lists), where protein target has a unique identifier of a UniProtKB accession, annotations also include protein family, aminmo acid sequences, and species.
* Known binding affinities (IC50, list, same order as in the protein target, could also be a list of lists, e.g. [[target 1, affinity], [target 2, affinity]])

### Nice to have

* 2d mol file (string/blob)
* known adverse effects (list)
* Other measured attributes, e.g. Kd (list)
* PAINS annotation(s) (boolean) -- does the molecule contain a functional group known to be "super active", i.e. shows up in assays even if the molecule itself is inactive
* protein target 3D structure(s) (list, same order and PT list)
* protein target associated disease(s) (list)

## Template

Include a link to your template file, or definition of the schema, e.g.:

```json
{"InChI-Key": {
    "pref_name": "Compound Name", "inchi_key": "InChI-Key",
    "synonyms": [
        {"synonym": "Synonym1", "source": "Source of Synonym1", "type": "Synonym Type"}],
    "structure": "Compound Structure",
    "properties": {"molecular_weight": "Molecular Weight", "formula": "Chemical Formula", "smiles": "SMILES Notation", "inchi": "InChI Notation", "logP": "number", "solubility": "number", "pKa": "number", "tpsa": "number", "h_bond_donors": "number", "h_bond_acceptors": "number", "rotatable_bonds": "number"},
    "targets": [
        {"target_name": "Target Name (e.g. Kinase Name)", "uniprot_accession": "UniProt Accession ID", "enzyme_family": "Enzyme Family", "activity_type": "Activity Type (e.g. IC50, binding affinity, etc)", "activity_value": "Activity Value", "unit": "Unit of Activity", "source": "Source of the activity, e.g. ChEMBL", "source_id": "identifier for the activity, e.g. paper DOI, or ChEMBL ID"}], 
    "references": [
        {"title": "string", "authors": "string", "journal": "string", "year": "number", "doi": "string"}]}}
```

## Possible Other Research Questions

### Which of my compounds have adverse effects?

The kinds of questions you might also consider here are:
* Can I predict if a new compound will have similar issues? 
* What are the most common adverse effects for my list of compounds?

Adverse effects, also known as side effects, are unintended and often harmful outcomes that occur when a drug or compound interacts with the body. These effects can range from mild, such as headaches or nausea, to severe, like organ damage or life-threatening allergic reactions. They also range in likelihood, from common to rare (or unknown), which reflect the number of people that have reported on these effects. Adverse effects can broken down into the following types:
 * Acute Effects: These occur shortly after exposure to a compound and are often reversible once the compound is removed.
 * Chronic Effects: These develop slowly over time with prolonged exposure and may be irreversible.
 * Idiosyncratic Reactions: Unpredictable reactions that occur in a small percentage of the population, often due to genetic differences.

Here are some suggestions for where to gather information on adverse effects:
* Many drug databases also have information on adverse effects (including ChEMBL), but specific databases also exist, such as SIDER (http://sideeffects.embl.de/)
* FDA Adverse Event Reporting System (FAERS): A public database that contains information on adverse event and medication error reports submitted to the FDA.
* European Medicines Agency (EMA) EudraVigilance: A system for managing and analyzing information on suspected adverse reactions to medicines authorized in the European Economic Area (EEA).
* ClinicalTrials.gov: A database of privately and publicly funded clinical studies conducted around the world. It includes information about adverse effects observed during clinical trials.
* WHO International Clinical Trials Registry Platform (ICTRP): Aims to ensure that a complete view of research is accessible to all those involved in health care decision making.
* Primary scientific literature, e.g. PubMed, EuroPMC, or Google Scholar.
* National Pharmacovigilance Centers: Many countries have centers dedicated to the collection and analysis of data on adverse drug reactions (ADRs).
* Yellow Card Scheme (UK): A system for collecting and monitoring information on suspected safety concerns or incidents involving medicines and medical devices.
* Electronic Health Records (EHRs): EHRs can be a rich source of real-world data on adverse effects, although accessing and standardizing this data can be challenging.
* Patient Registries: Disease-specific or treatment-specific registries collect long-term data on patients, including adverse effects.

### How soluble is my compound? 

The solubility of compounds plays a pivotal role in drug development, formulation, and overall effectiveness. Accurately predicting compound solubility can accelerate the development process, reduce costs, and improve therapeutic efficacy. To achieve robust predictive models, it is imperative to curate a comprehensive and high-quality dataset. This introduction will guide you through the types of data required and potential sources for compiling such a dataset.

To predict solubility, one must gather data on various chemical properties of the compounds. These properties include:
* Molecular Weight: The mass of a molecule, which can influence its solubility.
* LogP (Partition Coefficient): A measure of lipophilicity that affects solubility in different solvents.
* Hydrogen Bond Donors and Acceptors: The ability to form hydrogen bonds, which impacts solubility in polar solvents.
* Topological Polar Surface Area (TPSA): The surface area of polar atoms, which correlates with solubility in aqueous environments.
* pKa: The dissociation constant, which affects solubility at different pH levels.

It is also useful to collate experimental results, such as: 
* Solubility in Water: Baseline solubility in an aqueous environment.
* Solubility in Organic Solvents: Solubility in solvents such as ethanol, methanol, and DMSO.
* Solubility at Different pH Levels: Data indicating how solubility changes across the pH spectrum.
* Temperature-Dependent Solubility: Information on how solubility varies with temperature.
* Crystal structures: Information on crystalling forms as this can affect solubility.

There are lots of data sources out there, including public databases such as:
* ChemSpider: A free chemical structure database providing access to over 60 million structures, properties, and associated information.
* PubChem: A comprehensive resource for chemical properties, structures, and bioactivity data.
* DrugBank: A database containing detailed chemical, pharmacological, and pharmaceutical data on drugs and drug-like compounds.

You should also consider the scientific literature and even explore other predictive tools, such as SwissADME (swissadme.ch)

### Can I synthesise my compound utilising enzymes or low environmental impact methods?

The increasing awareness of environmental sustainability has led researchers to explore methods that reduce the ecological footprint of chemical syntheses. One pertinent research question is whether it is possible to synthesize compounds using enzymes or low environmental impact methods. 

To effectively synthesize compounds utilizing enzymes or low environmental impact methods, it is essential to gather a diverse range of data from reliable sources. Possible data sources include:
- Scientific Literature: Peer-reviewed journals and articles that provide detailed information on enzyme functions, reaction conditions, and case studies of green synthesis methods.
- Databases: Comprehensive databases such as the Protein Data Bank (PDB) for enzyme structures, Reaxys for chemical reactions, and PubChem for chemical properties and bioactivity data. Consider the overall chemistry being performed (e.g. the Enzyme Commission Number, Rhea and UniProtKB are good sources for this), and maybe even the reaction mechanism (e.g. M-CSA). You will also find kinetic data in databases such as BRENDA or SABIO-RK
- Patents: Patent databases can offer insights into novel synthesis methods and proprietary enzyme technologies.
- Environmental Impact Reports: Assessments and reports that evaluate the environmental outcomes of various synthesis methods, including energy consumption, waste generation, and overall sustainability.
- Experimental Data: Data generated from lab experiments, including reaction yields, enzyme activities, and conditions optimizing green synthesis pathways.
- Industry Reports: Documents from chemical industries that provide practical insights into large-scale applications of green synthesis methods.
- Regulatory Databases: Sources like the European Chemicals Agency (ECHA) that offer guidelines and compliance information relevant to environmentally friendly synthesis processes.

The kinds of data to gather could include:
- Enzyme Kinetics and Specificity: Information on enzyme activities, substrate specificity, and optimal reaction conditions.
- Chemical Properties: Data on melting points, boiling points, solubility, and stability of compounds.
- Synthesis Pathways: Detailed protocols and step-by-step procedures for synthesizing compounds using enzymes or other green methods.
- Environmental Metrics: Measurements of carbon footprint, waste by-products, energy efficiency, and overall environmental impact of the synthesis processes.
- Safety and Toxicity Data: Information on the potential hazards and safety measures required for handling chemicals and enzymes used in the synthesis.
- Economic Data: Cost analysis of raw materials, enzymes, and overall production expenses for green synthesis methods.

By gathering and curating data from these sources, researchers can better understand and implement sustainable synthesis methods, ultimately advancing the goal of a greener chemical industry.
