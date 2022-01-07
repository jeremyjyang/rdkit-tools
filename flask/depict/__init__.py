# Based on https://flask.palletsprojects.com/en/2.0.x/tutorial/factory/

import os,io,base64,logging

from flask import Flask,render_template
from markupsafe import escape

import rdkit.Chem
import rdkit.Chem.Draw
import rdkit.Chem.AllChem
#import collections
#from rdkit.Chem.Draw.MolDrawing import DrawingOptions

from rdktools import depict

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

def create_app(test_config=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        DATABASE=os.path.join(app.instance_path, 'flaskr.sqlite'),
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    # simple hello
    #@app.route('/hello')
    #def hello():
    #    return 'Hello, World!'

    # better hello
    @app.route('/hello/')
    @app.route('/hello/<name>')
    def hello(name=None):
      return render_template('hello.html', name=name)

    @app.route('/')
    def home():
        return 'Depict app (under construction).'

    @app.route('/mol2img/<smiles>', methods=['GET'])
    def mol2img(smiles):
      # show depiction for smiles
      mol = rdkit.Chem.MolFromSmiles(smiles)
      #mol = rdkit.Chem.AddHs(mol)
      rdkit.Chem.AllChem.Compute2DCoords(mol, clearConfs=True)
      # PIL.Image.Image
      img = depict.Utils.Mol2Image(mol, width=300, height=200, kekulize=True, wedgeBonds=True)
      img_bytearray = io.BytesIO()
      #img_bytearray.write(b'Content-type: image/png\n\n')
      img.save(img_bytearray, format='PNG')
      imgdata_b64 = base64.b64encode(img_bytearray.getvalue())
      logging.debug(f"SMILES: {escape(smiles)}")
      logging.debug(f"len(imgdata_b64): {len(imgdata_b64)}")
      return render_template('image.html', result=imgdata_b64)


    return app
