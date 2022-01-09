# Ref: https://flask.palletsprojects.com/en/2.0.x/tutorial/factory/

import os,io,base64,logging

from flask import Flask,render_template,request
#from flask_sqlalchemy import SQLAlchemy

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

    app.config.from_object(os.environ['APP_SETTINGS'])
    #app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
    #db = SQLAlchemy(app)

    #app.config.from_mapping(
    #    SECRET_KEY='dev',
    #    DATABASE=os.path.join(app.instance_path, 'flaskr.sqlite'),
    #)

    # PostgreSql db
    from . import db
    db.init_app(app)

    if test_config is None:
      # load the instance config, if it exists, when not testing
      app.config.from_pyfile('config.py', silent=True)
    else:
      # load the test config if passed in
      app.config.from_mapping(test_config)

    # ensure the instance folder exists
    if not os.path.exists(app.instance_path):
      try:
        os.makedirs(app.instance_path)
      except OSError:
        logging.error(f"Failed to create instance_path: {app.instance_path}")

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
      logging.info(f"SMILES: {escape(smiles)}")
      logging.debug(f"request.method: {request.method}")
      logging.debug(f"request.args.get('kekulize'): {request.args.get('kekulize')}")
      logging.debug(f"request.args.get('width'): {request.args.get('width')}")
      logging.debug(f"request.args.get('height'): {request.args.get('height')}")
      # show depiction for smiles
      mol = rdkit.Chem.MolFromSmiles(smiles)
      rdkit.Chem.AllChem.Compute2DCoords(mol, clearConfs=True)
      width = int(request.args.get('width')) if request.args.get('width') else 400
      height = int(request.args.get('height')) if request.args.get('height') else 400
      kekulize = bool(request.args.get('kekulize'))
      #img = depict.Utils.Mol2Image(mol, width=width, height=height, kekulize=kekulize, wedgeBonds=True)
      img = rdkit.Chem.Draw.MolToImage(mol, size=(width, height), kekulize=kekulize, wedgeBonds=True, fitImage=True) # PIL.Image.Image
      img_bytearray = io.BytesIO()
      img.save(img_bytearray, format='PNG')
      tmpfile = f"{app.instance_path}/tmp.png"
      fout = open(tmpfile, "wb+")
      fout.write(img_bytearray.getvalue())
      fout.close()
      imgdata_b64 = base64.b64encode(img_bytearray.getvalue()).decode("ascii")
      logging.debug(f"SEE: {tmpfile}")
      logging.debug(f"len(imgdata_b64): {len(imgdata_b64)}")
      #logging.debug(f"imgdata_b64: {imgdata_b64}")
      return render_template('image.html', result=imgdata_b64)


    return app
