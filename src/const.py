###########################
# Hard-coded inputs
###########################

levels = ["county", "tract"]

POP_CODE = "P0010001"

names = {"county": "counties", "tract": "tracts"}

state_codes = {
    'WA': '53', 'DE': '10', 'WI': '55', 'WV': '54', 'HI': '15',
    'FL': '12', 'WY': '56', 'NJ': '34', 'NM': '35', 'TX': '48',
    'LA': '22', 'NC': '37', 'ND': '38', 'NE': '31', 'TN': '47', 'NY': '36',
    'PA': '42', 'AK': '02', 'NV': '32', 'NH': '33', 'VA': '51', 'CO': '08',
    'CA': '06', 'AL': '01', 'AR': '05', 'VT': '50', 'IL': '17', 'GA': '13',
    'IN': '18', 'IA': '19', 'MA': '25', 'AZ': '04', 'ID': '16', 'CT': '09',
    'ME': '23', 'MD': '24', 'OK': '40', 'OH': '39', 'UT': '49', 'MO': '29',
    'MN': '27', 'MI': '26', 'RI': '44', 'KS': '20', 'MT': '30', 'MS': '28',
    'SC': '45', 'KY': '21', 'OR': '41', 'SD': '46'
}

congressional_districts = {
    'WA': 10, 'DE': 1, 'WI': 8, 'WV': 3, 'HI': 2,
    'FL': 27, 'WY': 1, 'NJ': 12, 'NM': 3, 'TX': 36,
    'LA': 6, 'NC': 13, 'ND': 1, 'NE': 3, 'TN': 9, 'NY': 27,
    'PA': 18, 'AK': 1, 'NV': 4, 'NH': 2, 'VA': 11, 'CO': 7,
    'CA': 53, 'AL': 7, 'AR': 4, 'VT': 1, 'IL': 18, 'GA': 14,
    'IN': 9, 'IA': 4, 'MA': 9, 'AZ': 9, 'ID': 2, 'CT': 5,
    'ME': 2, 'MD': 8, 'OK': 5, 'OH': 16, 'UT': 4, 'MO': 8,
    'MN': 8, 'MI': 14, 'RI': 2, 'KS': 4, 'MT': 1, 'MS': 4,
    'SC': 7, 'KY': 6, 'OR': 5, 'SD': 1
}

skips = {
    ('WA','tract'), ('WA','county'), ('DE','tract'), ('DE','county'), ('WI','tract'),
    ('WI','county'), ('HI','tract'), ('HI','county'), ('FL','tract'), ('FL','county'),
    ('WY','tract'), ('WY','county'), ('NJ','tract'), ('NJ','county'), ('TX','tract'),
    ('TX','county'), ('LA','tract'), ('LA','county'), ('NC','tract'), ('NC','county'),
    ('ND','tract'), ('ND','county'), ('TN','tract'), ('TN','county'), ('NY','tract'),
    ('NY','county'), ('PA','tract'), ('PA','county'), ('AK','tract'), ('AK','county'),
    ('NV','county'), ('NH','county'), ('VA','tract'), ('VA','county'), ('CO','tract'),
    ('CO','county'), ('CA','tract'), ('CA','county'), ('AL','tract'), ('VT','tract'),
    ('VT','county'), ('IL','tract'), ('IL','county'), ('GA','tract'), ('GA','county'),
    ('IN','tract'), ('IN','county'), ('IA','tract'), ('MA','tract'), ('MA','county'),
    ('AZ','tract'), ('AZ','county'), ('CT','tract'), ('CT','county'), ('ME','county'),
    ('MD','tract'), ('MD','county'), ('OK','tract'), ('OH','tract'), ('OH','county'),
    ('UT','county'), ('MO','tract'), ('MO','county'), ('MN','tract'), ('MN','county'),
    ('MI','tract'), ('MI','county'), ('RI','tract'), ('RI','county'), ('KS','tract'),
    ('MT','tract'), ('MT','county'), ('SC','tract'), ('SC','county'), ('KY','tract'),
    ('KY','county'), ('OR','tract'), ('OR','county'), ('SD','tract'), ('SD','county')
}