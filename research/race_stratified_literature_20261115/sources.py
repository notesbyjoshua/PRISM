import json
import re
import unicodedata
import xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path

import pandas as pd


def normalize(text):
    text = unicodedata.normalize('NFKD', str(text)).encode('ascii', 'ignore').decode()
    return re.sub(r'[^a-z0-9]', '', text.lower())


def ontology(path):
    terms, current = {}, None
    for line in Path(path).read_text().splitlines():
        if line == '[Term]':
            current = {'parents': [], 'synonyms': [], 'alt_ids': []}
        elif line.startswith('['):
            current = None
        elif current is not None:
            if line.startswith('id: HP:'):
                current['id'] = line[4:]
                terms[current['id']] = current
            elif line.startswith('name: '):
                current['name'] = line[6:]
            elif line.startswith('def: '):
                current['definition'] = line[5:]
            elif line.startswith('is_a: HP:'):
                current['parents'].append(line.split()[1])
            elif line.startswith('synonym: '):
                match = re.match(r'synonym: "(.*?)"', line)
                if match:
                    current['synonyms'].append(match.group(1))
            elif line.startswith('alt_id: '):
                current['alt_ids'].append(line[8:])
            elif line == 'is_obsolete: true':
                current['obsolete'] = True
    aliases = {alt: key for key, term in terms.items() for alt in term['alt_ids']}
    names = {}
    for key, term in terms.items():
        if term.get('obsolete'):
            continue
        for value in [term.get('name', '')] + term['synonyms']:
            names.setdefault(normalize(value), set()).add(key)

    @lru_cache(maxsize=None)
    def ancestors(key):
        key = aliases.get(key, key)
        return frozenset([key]).union(*(ancestors(p) for p in terms.get(key, {}).get('parents', [])))

    return terms, names, ancestors


def annotation_tables(directory):
    hpo = pd.read_csv(directory / 'phenotype.hpoa', comment='#', sep='\t', dtype=str).fillna('')
    hpo = hpo.rename(columns={'database_id': 'disease_id', 'hpo_id': 'term_id'})
    hpo['source'] = 'HPO_annotations_2026-09-02'
    hpo['negated'] = hpo.qualifier.str.contains('NOT')
    root = ET.parse(directory / 'orphanet_hpo.xml').getroot()
    rows = []
    for disease in root.findall('.//HPODisorderSetStatus/Disorder'):
        for association in disease.findall('./HPODisorderAssociationList/HPODisorderAssociation'):
            frequency = association.findtext('./HPOFrequency/Name', '')
            rows.append({'disease_id': 'ORPHA:' + disease.findtext('OrphaCode'), 'disease_name': disease.findtext('Name'),
                         'term_id': association.findtext('./HPO/HPOId'), 'reference': '', 'evidence': 'Orphanet expert annotation',
                         'frequency': frequency, 'onset': '', 'sex': '', 'modifier': '',
                         'source': 'Orphadata_HPO_2026-06-23', 'negated': frequency.startswith('Excluded')})
    return pd.concat([hpo, pd.DataFrame(rows)], ignore_index=True)


def nomenclature(path):
    root = ET.parse(path).getroot()
    return [{'orpha': 'ORPHA:' + item.findtext('OrphaCode'), 'name': item.findtext('Name'),
             'synonyms': [node.text for node in item.findall('./SynonymList/Synonym') if node.text],
             'type': item.findtext('./DisorderType/Name', ''),
             'omim': [ref.findtext('Reference') for ref in item.findall('./ExternalReferenceList/ExternalReference')
                      if ref.findtext('Source') == 'OMIM']} for item in root.findall('./DisorderList/Disorder')]
