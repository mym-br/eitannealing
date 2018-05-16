#include <QCommandLineParser>
#include "parametersparser.h"
#include "problem.h"

CommandLineParseResult parseCommandLine(QCommandLineParser &parser, EitAnnealingArgs *params, QString *errorMessage)
{
    parser.setSingleDashWordOptionMode(QCommandLineParser::ParseAsLongOptions);

    parser.addPositionalArgument("mesh","Input mesh file.");
	parser.addPositionalArgument("currents", "Input injected currents file.");
	parser.addPositionalArgument("tensions", "Input tension measurements file.");
	parser.addPositionalArgument("regularisation", "Weight for gradient norm regularisation.");
	parser.addPositionalArgument("electrodevar", "Weight for electrode variance penalization.");
	parser.addPositionalArgument("kt", "Initial temperature.");
	
	const QCommandLineOption gmeshOpt("gmesh", "Address for gmesh socket communication.", "address"); parser.addOption(gmeshOpt);
	const QCommandLineOption seedOpt("seed", "Seed for the random number generator.", "number"); parser.addOption(seedOpt);
	const QCommandLineOption outputOpt("output", "Output mesh file name.", "number"); parser.addOption(outputOpt);
	const QCommandLineOption peparamOpt("peparam", "Partial cost evaluation probability parameter.", "probability"); parser.addOption(peparamOpt);
	const QCommandLineOption gndOpt("ground", "Ground node.", "number"); parser.addOption(gndOpt);
	const QCommandLineOption currentsOutOpt("currentsout", "Input removed currents file.", "number"); parser.addOption(currentsOutOpt);
	const QCommandLineOption calibrationOpt("calibration", "Calibration mode: 1 - single coefficient for electrodes; 2 - individual coefficients for electrodes.", "number"); parser.addOption(calibrationOpt);
	const QCommandLineOption ktOpt("kt", "Initial temperature.", "number"); parser.addOption(ktOpt);
	const QCommandLineOption regularizationOpt("regularization", "Weigth for the gradirnt norm regularization.", "number"); parser.addOption(regularizationOpt);
	const QCommandLineOption electrodeVarOpt("elecvar", "Weigth for the electrode variance regularization.", "number"); parser.addOption(electrodeVarOpt);

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
	if (parser.isSet(currentsOutOpt)) params->inputCurrentsOut = parser.value(currentsOutOpt);
	if (parser.isSet(gmeshOpt)) params->gmeshAddress = parser.value(gmeshOpt); else params->gmeshAddress = "127.0.0.1:44202";
	if (parser.isSet(outputOpt)) params->outputMesh = parser.value(outputOpt); else params->outputMesh = "solution.msh";
	if (parser.isSet(seedOpt)) params->setSeed(parser.value(seedOpt).toULong());
	if (parser.isSet(peparamOpt)) params->peParam = parser.value(peparamOpt).toFloat(); else params->peParam = 0.875f;
	if (parser.isSet(gndOpt)) params->ground = parser.value(gndOpt).toInt(); else params->ground = -1;
	if (parser.isSet(calibrationOpt)) params->calibrationMode = parser.value(calibrationOpt).toInt();
	if (parser.isSet(ktOpt)) params->kt = parser.value(ktOpt).toDouble(); else params->kt = 0.05;
	if (parser.isSet(regularizationOpt)) params->regularizationFactor = parser.value(regularizationOpt).toDouble(); else params->regularizationFactor = 30;
	if (parser.isSet(electrodeVarOpt)) params->electrodevar = parser.value(electrodeVarOpt).toDouble(); else params->electrodevar = 0;

    return CommandLineOk;
}