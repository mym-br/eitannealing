#include <QCommandLineParser>
#include "parametersparser.h"

CommandLineParseResult parseCommandLine(QCommandLineParser &parser, EitAnnealingArgs *params, QString *errorMessage)
{
    parser.setSingleDashWordOptionMode(QCommandLineParser::ParseAsLongOptions);

    parser.addPositionalArgument("mesh","Input mesh file.");
	parser.addPositionalArgument("currents", "Input injected currents file.");
	parser.addPositionalArgument("tensions", "Input tension measurements file.");
	
	const QCommandLineOption gmeshOpt("gmesh", "Address for gmesh socket communication.", "address"); parser.addOption(gmeshOpt);
	const QCommandLineOption seedOpt("seed", "Seed for the random number generator.", "number"); parser.addOption(seedOpt);

    const QCommandLineOption helpOption = parser.addHelpOption();
    const QCommandLineOption versionOption = parser.addVersionOption();

    if (!parser.parse(QCoreApplication::arguments())) {
        *errorMessage = parser.errorText();
        return CommandLineError;
    }

    if (parser.isSet(versionOption))
        return CommandLineVersionRequested;

    if (parser.isSet(helpOption))
        return CommandLineHelpRequested;

    const QStringList positionalArguments = parser.positionalArguments();
    if (positionalArguments.isEmpty()) {
        *errorMessage = "No argument specified.";
        return CommandLineError;
    }
	if (positionalArguments.size() < 3) {
		*errorMessage = "Too few arguments.";
		return CommandLineError;
	}
    if (positionalArguments.size() > 3) {
        *errorMessage = "Too many arguments specified.";
        return CommandLineError;
    }
    params->inputMesh = positionalArguments.at(0);
	params->inputCurrents = positionalArguments.at(1);
	params->inputTensions = positionalArguments.at(2);

	if (parser.isSet(gmeshOpt)) params->gmeshAddress = parser.value(gmeshOpt);
	else params->gmeshAddress = "127.0.0.1:44202";
	if (parser.isSet(seedOpt)) params->setSeed(parser.value(seedOpt).toULong());

    return CommandLineOk;
}