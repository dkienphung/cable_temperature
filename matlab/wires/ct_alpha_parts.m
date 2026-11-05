function parts = ct_alpha_parts()
%CT_ALPHA_PARTS  The 337 Alpha Wire records of data/alphawire_hookup.json (cell array of structs).
db = ct_json('alphawire_hookup.json');
parts = ct_as_cell(db.parts);
end
