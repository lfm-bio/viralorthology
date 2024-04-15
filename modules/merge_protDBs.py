import os

def merge(db):
    new_db = open(db, 'a', encoding='utf-8')
    first_round_db_name = f'first-round-{db}'
    try:
        first_round_db = open(first_round_db_name, encoding='utf-8')
    except:
        return False #there werent filtered genomes
    for line in first_round_db:
        new_db.write(line)
    new_db.close()
    first_round_db.close()
    os.remove(first_round_db_name)

def main():
    merge('protDB.db')
    merge('protDB_OF.db')