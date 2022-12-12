from argparse import ArgumentParser
from coreapi import Client
import os

def get_input():
    parser = ArgumentParser()
    parser.add_argument(
        "tf", type = str, help="TF for which you want to query motifs from JASPAR"
    )
    return parser.parse_args()

if __name__ == "__main__":
    client = Client()
    document = client.get('https://jaspar.genereg.net/api/v1/docs/')

    tf = get_input().tf
    data = client.action(document, ['matrix','list'], params = {
        'search': tf,
        'tax_id': '9606',
        'release': 'JASPAR2022'
    })

    if data['count'] != 0:
        os.makedirs(tf)
        for motif in data['results']:
            if tf in motif['name']:
                id = motif['matrix_id']
                query = str(f"https://jaspar.genereg.net/api/v1/matrix/{id}/?format=jaspar")
                pfm = client.get(query)
                with open(f"{tf}/{id}.jaspar", 'w') as outfile:
                    outfile.writelines(pfm)
    else:
        print(f"No motifs in JASPAR for {tf}")
