import psycopg2,psycopg2.extras
import click
from flask import current_app, g
from flask.cli import with_appcontext


def get_db():
    if 'db' not in g:
      #dsn = (f"host='{dbhost}' port='{dbport}' dbname='{dbname}' user='{dbusr}' password='{dbpw}'")
      #dbcon = psycopg2.connect(dsn)
      g.db = sqlite3.connect(current_app.config['DATABASE'])
      dbcon.cursor_factory = psycopg2.extras.DictCursor
    return g.db


def close_db(e=None):
    db = g.pop('db', None)
    if db is not None:
      db.close()

def open_db():
    db = get_db()
    with current_app.open_resource('describe.sql') as f:
      db.executescript(f.read().decode('utf8'))

@click.command('open-db')
@with_appcontext
def open_db_command():
    """Open the existing db and describe."""
    open_db()
    click.echo('Opened the database.')

def init_app(app):
    app.teardown_appcontext(close_db)
    app.cli.add_command(open_db_command)
