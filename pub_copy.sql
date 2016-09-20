CREATE TABLE "photometry" (
	`id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	`source_id`	INTEGER NOT NULL,
	`band`	TEXT NOT NULL,
	`magnitude`	REAL NOT NULL,
	`magnitude_unc`	REAL,
	`system`	INTEGER,
	`telescope_id`	INTEGER,
	`instrument_id`	INTEGER,
	`publication_id`	INTEGER,
	`epoch`	TEXT,
	`comments`	TEXT
);
CREATE TABLE "systems" (
	`id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	`name`	TEXT NOT NULL
);
CREATE TABLE "radial_velocities" (
	`id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	`source_id`	INTEGER NOT NULL,
	`radial_velocity`	REAL NOT NULL,
	`radial_velocity_unc`	REAL,
	`spectrum_id`	INTEGER,
	`publication_id`	INTEGER
);
CREATE TABLE "modes" (
	`id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	`mode`	TEXT NOT NULL
);
CREATE TABLE "instruments" (
	`id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	`name`	TEXT NOT NULL,
	`publication_id`	INTEGER
);
CREATE TABLE "telescopes" (
	`id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	`name`	TEXT NOT NULL,
	`publication_id`	INTEGER
);
CREATE TABLE "publications" (
	`id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	`bibtex`	TEXT,
	`shortname`	TEXT,
	`DOI`	TEXT,
	`description`	TEXT
);
CREATE TABLE "spectra" (
	`id`	INTEGER,
	`source_id`	INTEGER,
	`spectrum`	SPECTRUM,
	`wavelength_units`	TEXT,
	`flux_units`	TEXT,
	`wavelength_order`	INTEGER,
	`regime`	TEXT,
	`publication_id`	INTEGER,
	`obs_date`	TEXT,
	`instrument_id`	INTEGER,
	`telescope_id`	INTEGER,
	`mode_id`	INTEGER,
	`filename`	TEXT,
	`comments`	TEXT,
	PRIMARY KEY(id)
);
CREATE TABLE "parallaxes" (
	`id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	`source_id`	INTEGER NOT NULL,
	`parallax`	REAL NOT NULL,
	`parallax_unc`	REAL,
	`publication_id`	INTEGER,
	`adopted`	INTEGER,
	`comments`	TEXT
);
CREATE TABLE "proper_motions" (
	`id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	`source_id`	INTEGER NOT NULL,
	`proper_motion_ra`	REAL NOT NULL,
	`proper_motion_ra_unc`	REAL,
	`proper_motion_dec`	REAL NOT NULL,
	`proper_motion_dec_unc`	REAL,
	`spectrum_id`	INTEGER,
	`publication_id`	INTEGER,
	`comments`	TEXT
);
CREATE TABLE "spectral_types" (
	`id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	`source_id`	INTEGER NOT NULL,
	`spectral_type`	REAL NOT NULL,
	`spectral_type_unc`	REAL,
	`gravity`	TEXT,
	`suffix`	TEXT,
	`regime`	TEXT,
	`publication_id`	INTEGER,
	`comments`	TEXT,
	`adopted`	BOOLEAN
);
CREATE TABLE "sources" (
	`id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	`ra`	REAL,
	`dec`	REAL,
	`designation`	TEXT,
	`publication_id`	INTEGER,
	`comments`	TEXT,
	`unum`	TEXT,
	`shortname`	TEXT,
	`names`	TEXT,
	`components`	TEXT,
	`companions`	TEXT
);
CREATE TABLE ignore (id INTEGER PRIMARY KEY, id1 INTEGER, id2 INTEGER, tablename TEXT);
