# This model is adapted from previous work from Katarzyna Wardega and Wendy Mendoza

from keras import models, layers
from tensorflow.keras.utils import plot_model
from sklearn.metrics import (
        accuracy_score,
        classification_report,
        ConfusionMatrixDisplay,
        confusion_matrix,
        RocCurveDisplay,
)


from config import Configuration


CUTOUT_WIDTH = Configuration.CUTOUT_WIDTH

def train():
    # create model
    model = models.Sequential()

    # define input shape
    model.add(layers.Input(shape=(CUTOUT_WIDTH, CUTOUT_WIDTH, 2)))

    # add model layers
    model.add(layers.Conv2D(10, kernel_size=3, activation='relu'))
    model.add(layers.Conv2D(5, kernel_size=3, activation='relu'))
    model.add(layers.MaxPooling2D(pool_size=(3,3)))
    model.add(layers.Dropout(0.25))
    model.add(layers.Conv2D(3, kernel_size=3, activation='relu'))
    model.add(layers.MaxPooling2D(pool_size=(2,2)))
    model.add(layers.Flatten())
    model.add(layers.Dense(10, activation='relu'))
    model.add(layers.Dropout(0.5))
    model.add(layers.Dense(50, activation='relu'))
    model.add(layers.Dropout(0.3))
    model.add(layers.Dense(2, activation='softmax'))

    # show summary
    model.summary()

    # save summary as image
    plot_model(model, to_file='model_summary.png', show_shapes=True, show_layer_names=True)

    # Set up training
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

    # Train
    history = model.fit(
            data,
            labels,
            validation_data=(data_test, labels_test), 
            epochs=30)

    # Save model
    model.save("toros_cnn_model.h5")

    # plot metrics
    # plot accuracy vs epoch
    plt.plot(history.history['accuracy'])
    plt.plot(history.history['val_accuracy'])
    plt.title('model accuracy')
    plt.ylabel('accuracy')
    plt.xlabel('epoch')
    plt.legend(['train', 'test'], loc='upper left')
    plt.savefig('accuracy_seq.png')
    plt.show()

    report = classification_report(y_test, y_pred, target_names=['bogus', 'real'],
                                   output_dict=True)
    print(classification_report(y_test, y_pred, target_names=['bogus', 'real']))

    # Confusion Matrix
    cm = confusion_matrix(y_test, y_pred)
    cm_display = ConfusionMatrixDisplay(cm).plot()
    plt.title(f"Confusion Matrix of Test Set\n bogus: {report['bogus']['support']} real: {report['real']['support']}\n Accuracy Score: {accuracy_score(y_test, y_pred):.4f} ")
    plt.savefig("ConfusionMatrix.png", dpi=300, bbox_inches="tight")
    plt.show()
