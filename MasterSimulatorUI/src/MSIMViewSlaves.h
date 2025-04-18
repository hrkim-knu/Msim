#ifndef MSIMVIEWSLAVES_H
#define MSIMVIEWSLAVES_H

#include <QWidget>

#include <MSIM_Project.h>

namespace QPW {
	class VariantPropertyManager;
}

namespace Ui {
	class MSIMViewSlaves;
}

namespace IBK {
	class Path;
}

namespace MASTER_SIM {
	class ModelDescription;
}

namespace BLOCKMOD {
	class BlockItem;
}

class MSIMBlockEditorDialog;
class QTableWidgetItem;
class MSIMExportConnectionGraphDialog;

/*! The view containing FMU and slave definition tables. */
class MSIMViewSlaves : public QWidget {
	Q_OBJECT

public:

	enum SelectionState {
		SS_SlaveSelected,
		SS_ConnectorSelected,
		SS_NothingSelected
	};


	/*! C'tor */
	explicit MSIMViewSlaves(QWidget *parent = nullptr);
	/*! D'tor */
	~MSIMViewSlaves();

	/*! Looks up the block item corresponding to the given slave and opens the block editor. */
	void editBlockItem(const QString & slaveName);

	/*! Opens print/export dialog and exports/prints scene. */
	void printScene();

	/*! Extract the FMU with the given absolute file path, and determined properties for analysis.
		If successful, returns true and the modelDesc data structure can be accessed to provide the
		required information.
		The msgLog is a log of the operations done and can be shown in a plain text edit field.
	*/
	static bool extractFMUAndParseModelDesc(const IBK::Path & fmuFilePath,
									QString & msgLog,
									MASTER_SIM::ModelDescription & modelDesc,
									QPixmap & modelPixmap);

public slots:
	/*! Connected to MSIMProjectHandler::modified() */
	void onModified(unsigned int modificationType, void * data );

signals:
	/*! Emitted when a new slave has been added.
		The second argument passed is the **absolute file path** to the slave (FMU or tsv/csv table).
	*/
	void newSlaveAdded(const QString & slaveName, const QString & fmuFilePath);

private slots:
	void on_toolButtonAddSlave_clicked();

	void on_toolButtonRemoveSlave_clicked();

	void on_tableWidgetSlaves_cellChanged(int row, int column);

	void on_checkBoxRelativeFMUPaths_toggled(bool checked);

	void on_tableWidgetSlaves_currentCellChanged(int , int , int , int );

	/*! Connected to the scene's block triggering action (usually double-click on block item). */
	void onBlockActionTriggered(const BLOCKMOD::BlockItem * blockItem);

	/*! Function is called from single-shot timer once a block-action trigger slot has finished.
		Creates an undo-action for a block change.
	*/
	void onBlockEditingCompleted();

	/*! Called from sceneManager() whenever user has finished creating a connection. */
	void onNewConnectionCreated();

	/*! Called from sceneManager() whenever a block or connector has been moved. */
	void onNetworkGeometryChanged();

	/*! Called from sceneManager() whenever a block has been selected. */
	void onBlockSelected(const QString & blockName);

	/*! Called from sceneManager() whenever a connector has been selected. */
	void onConnectorSelected(const QString & sourceSocketName, const QString & targetSocketName);

	/*! Called from sceneManager() whenever selection has changed and is now empty (there is nothing selected) */
	void onSelectionCleared();

	/*! User has edited a parameter of a slave. */
	void on_widgetProperties_itemChanged(QTableWidgetItem *item);

	void on_widgetConnectors_itemChanged(QTableWidgetItem *item);

	void on_pushButtonDeleteConnection_clicked();

	void on_doubleSpinBoxLinewidth_valueChanged(double arg1);

	void on_pushButtonSelectColor_colorChanged();

	void on_checkBoxShowEquations_clicked(bool checked);


	void on_textEditDescription_editingFinished();

	void on_pushButtonExportImage_clicked();

private:
	/*! Updates the table with all slaves defined for this simulation scenario. */
	void updateSlaveTable();

	/*! Populates the slave-specific parameter table which shows all parameters set by MasterSim.
		\param slaveIndex An index of a slave/simulator in the project. -1 means non selected/available.
	*/
	void updateSlaveParameterTable(unsigned int slaveIndex);

	/*! Populates the table widget with offset, scaleFactor values and updates the linewidth spinBox
		\param edge the current
	*/
	void updateGraphProperties();

	/*! Updates Ui property page. Shows slave properties, connector properties or nothing based on current selection state.
	 */
	void updatePropertyStackedWidget(SelectionState selectionState);

	Ui::MSIMViewSlaves				*m_ui;

	/*! The manager for the properties. */
	QPW::VariantPropertyManager		*m_variantManager = nullptr;

	MSIMExportConnectionGraphDialog	*m_exportConnectionGraphDialog = nullptr;

	MSIMBlockEditorDialog			*m_blockEditorDialog;

	int								m_selectedEdgeIdx = -1;
};

#endif // MSIMVIEWSLAVES_H
