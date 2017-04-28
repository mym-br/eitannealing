#ifndef PARAMETERSPARSER_H
#define PARAMETERSPARSER_H

class QCommandLineParser;
#include<QString>

struct EitAnnealingArgs {
	EitAnnealingArgs() : seedSpecified(false) {}

    QString inputMesh, inputCurrents, inputTensions;
	QString outputMesh;
	QString gmeshAddress;
	float peParam;
	int ground;

	bool isSeedSpecified() { return seedSpecified; }
	unsigned long getSeed() { return seed; }
	void setSeed(unsigned long _seed) { seed = _seed; seedSpecified = true; }

private:
	unsigned long seed;
	bool seedSpecified;
};

enum CommandLineParseResult
{
    CommandLineOk,
    CommandLineError,
    CommandLineVersionRequested,
    CommandLineHelpRequested
};

CommandLineParseResult parseCommandLine(QCommandLineParser &parser, EitAnnealingArgs *params, QString *errorMessage);
CommandLineParseResult parseOptionsFile(QString fileName, EitAnnealingArgs *params, QString *errorMessage);

#endif // PARAMETERSPARSER_H
