#ifndef MSIMMainWindow_H
#define MSIMMainWindow_H

#include <QMainWindow>
#include <QUndoStack>
#include <QProcess>
#include <QStackedWidget>

#include <map>

#include <MSIM_ModelDescription.h>

#include "MSIMProjectHandler.h"

namespace Ui {
	class MSIMMainWindow;
}

class MSIMWelcomeScreen;
class MSIMLogWidget;
class MSIMThreadBase;
class MSIMPreferencesDialog;
class MSIMViewSlaves;
class MSIMViewConnections;
class MSIMViewSimulation;
class MSIMAboutDialog;
class MSIMButtonBar;
class MSIMPostProcHandler;

/*! Main window class. */
class MSIMMainWindow : public QMainWindow {
	Q_OBJECT
public:

	/*! Returns a pointer to the MSIMMainWindow instance.
		Only access this function during the lifetime of the
		MSIMMainWindow instance.
	*/
	static MSIMMainWindow & instance();

	/*! Adds an undo command to the global undo stack.
		Ownership of the command object will be transferred to the stack.
	*/
	static void addUndoCommand(QUndoCommand * command);

	/*! Add a model description to the global list of known model descriptions in UI. */
	static void addModelDescription(const IBK::Path & fmuPath, const MASTER_SIM::ModelDescription & modelDesc, const QPixmap & modelPixmap);


	/*! Default MSIMMainWindow constructor. */
	explicit MSIMMainWindow(QWidget *parent = nullptr, Qt::WindowFlags flags = Qt::WindowFlags());

	/*! Default destructor. */
	~MSIMMainWindow();

	/*! Public access function to save project file (called from simulation view).
		\return Returns true if project was saved and project handler
				has now an unmodified project with valid project filename.
	*/
	bool saveProject();

	/*! Returns model descriptions from FMUs. */
	const std::map<IBK::Path, MASTER_SIM::ModelDescription>	& modelDescriptions() const { return m_modelDescriptions; }

	/*! Returns model description for a given slave name.
		Throws an IBK::Exception if either slave name is invalid, or FMU model description is not in the map.
	*/
	const MASTER_SIM::ModelDescription & modelDescription(const std::string & slaveName) const;
	const QPixmap & modelPixmap(const std::string & slaveName) const;

protected:
	/*! Checks if project file has been changed by external application. */
	void changeEvent(QEvent *event);
	/*! Does the confirm-saving-before-close stuff. */
	void closeEvent(QCloseEvent * event);

private slots:
	void on_actionFileNew_triggered();
	void on_actionFileOpen_triggered();
	void on_actionFileSave_triggered();
	void on_actionFileSaveAs_triggered();
	void on_actionFileReload_triggered();
	void on_actionFileExport_triggered();
	void on_actionFileClose_triggered();
	void on_actionFilePrint_triggered();
	void on_actionFileOpenProjectDir_triggered();
	void on_actionFileQuit_triggered();

	void on_actionEditTextEditProject_triggered();
		void on_actionEditParseFMUs_triggered();
	void on_actionEditPreferences_triggered();
	/*! Post-proc opening triggered. */
	void on_actionEditOpenPostProc_triggered();

	void on_actionHelpAboutQt_triggered();
	void on_actionHelpAboutMasterSim_triggered();
	void on_actionHelpMasterSimWebpage_triggered();
	void on_actionHelpOnlineManual_triggered();
	void on_actionHelpReport_triggered();
	void on_actionHelpLinuxDesktopIntegration_triggered();

	/*! Triggered when a recent file menu entry was clicked. */
	void onActionOpenRecentFile();

	/*! Triggered when a language menu entry was clicked. */
	void onActionSwitchLanguage();

	/*! Updates the state of all actions based on the current condition of the project.
		This slot is connected to the signal updateActions() from MSIMProjectHandler.
	*/
	void onUpdateActions();

	/*! Updates the menu entries in the 'Recent Projects' submenu.
		This is a slot because we need to update the menu with the actions
		a bit delayed. When the user clicks on a recent project menu entry, the
		loadProject() function is indirectly called which in turn calls
		updateRecentProjects(). Since the menu actions in the recent projects
		menu are deleted, this would mean that the action currently process is
		being deleted - causing a crash. Therefore we don't call updateRecentProjects()
		directly, but via a QTimer::singleShot() and thus ensure that the
		action handler function is completed before the action is touched.
	*/
	void onUpdateRecentProjects();

	/*! Opens a project with filename.
		Called from onActionOpenRecentFile() and from welcome screen.
	*/
	void onOpenProjectByFilename(const QString & filename);

	/*! Copies an example project (and its directory content) to a new directory (interactive) and
		opens the copied example project.
		Called only from welcome screen.
	*/
	void onOpenExampleByFilename(const QString & filename);

	/*! Slave view activated. */
	void on_actionViewSlaves_toggled(bool arg1);
	/*! Connvection view activated. */
	void on_actionViewConnections_toggled(bool arg1);
	/*! Simulation options view activated. */
	void on_actionViewSimulation_toggled(bool arg1);

	/*! Direct simulation start triggered (via shortcut). */
	void on_actionStartSimulation_triggered();

	/*! Called from slave view whenever an FMU has been added successfully (including successful
		parsing of the model description.
	*/
	void onNewSlaveAdded(const QString & slaveName, const QString & fullFMUPath);



private:
	/*! Sets up all dock widgets with definition lists. */
	void setupDockWidgets();

	/*! Updates the window title. */
	void updateWindowTitle();

	/*! Extracts all FMUs currently being referenced by the project and attempts to read their model description files.
	*/
	void extractFMUsAndParseModelDesc(QString & msg);

	/*! Exports the current project to the selected exportFilePath path.
		Opens an (internationalized) message box on the first error encountered.
		\return Returns true on success, false on error.
	*/
	bool exportProjectPackage(const QString & exportFilePath);

	/*! Imports the project package to the selected import directory.
		Opens an (internationalized) message box on the first error encountered.
		\param packageFilePath Full path to package file.
		\param targetDirectory Target dir to extract package content into (no subdirectory is created)
		\param projectFilePath If package contains a project file, this file name (full path) will be stored
		\param packageContainsMSIM True if package is an MSIM package and should contain an msim project file.
		\return Returns true on success, false on error.
	*/
	bool importProjectPackage(const QString & packageFilePath, const QString & targetDirectory,
							  QString & projectFilePath, bool packageContainsMSIM);

	/*! Creates a thumbnail-image of the current project sketch. */
	void saveThumbNail();

	/*! Adds another language setting action, when the corresponding language files exist. */
	void addLanguageAction(const QString & langId, const QString & actionCaption);

	/*! Helper function to remove a directory recursivly.
		\note in Qt5, see QDir::removeRecursively()
	*/
	static bool removeDirRecursively(const QString & dirName);

	/*! Global pointer to main window instance.
		Initialized in the constructor of MSIMMainWindow and
		reset to NULL in the destructor. So be sure that the main window
		exists before accessing MSIMMainWindow::instance()
	*/
	static MSIMMainWindow		*m_self;

	/*! Main user interface pointer. */
	Ui::MSIMMainWindow			*m_ui;
	/*! The global undo stack in the program. */
	QUndoStack					*m_undoStack;
	/*! Menu for the recent projects entries. */
	QMenu						*m_recentProjectsMenu;
	/*! List with action objects for each recent project in the main menu. */
	QList<QAction*>				m_recentProjectActions;
	/*! List of language actions shown in the languages menu. */
	QList<QAction*>				m_languageActions;
	/*! The project handler that manages the actual project. */
	MSIMProjectHandler			m_projectHandler;

	/*! Model description data for all FMUs currently parsed.
		This map is updated in on_actionEditParseFMUs_triggered().
		Key is the absolute file path to the FMU, value is the model description.
		The values are updated whenever the fmu is parsed.
	*/
	std::map<IBK::Path, MASTER_SIM::ModelDescription>	m_modelDescriptions;
	std::map<IBK::Path, QPixmap>						m_modelPixmaps;

	/*! The welcome screen. */
	MSIMWelcomeScreen			*m_welcomeScreen;

	/*! User preferences. */
	MSIMPreferencesDialog		*m_preferencesDialog;

	/*! Widget for logging content. */
	MSIMLogWidget				*m_logWidget;

	QStackedWidget				*m_stackedWidget;

	MSIMViewSlaves				*m_viewSlaves;
	MSIMViewConnections			*m_viewConnections;
	MSIMViewSimulation			*m_viewSimulation;

	MSIMAboutDialog				*m_aboutDialog;

	MSIMButtonBar				*m_buttonBar;

	/*! Handles spawning/activating of PostProc executable. */
	MSIMPostProcHandler			*m_postProcHandler;
};

#endif // MSIMMainWindow_H
