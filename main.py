from dbHandler import Handler
import os

json_path = os.path.join(os.path.dirname(__file__), "Sars-cov-2.json")

handler = Handler(json_path)

handler.display_statistics()
