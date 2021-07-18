#python3 t.py >AB981583.gb.feature.type.source.log
#cat AB981583.gb NC_005816.gb >combine.two.items.gb
python3 ../parse_gb_file.v2.py -gb combine.two.items.gb -of hehe
