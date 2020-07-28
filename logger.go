package main

import (
	"os"
	"time"

	"go.uber.org/zap"
	"go.uber.org/zap/zapcore"
	"gopkg.in/natefinch/lumberjack.v2"
)

// NewEncoderConfig as name says
func NewEncoderConfig() zapcore.EncoderConfig {
	return zapcore.EncoderConfig{
		// Keys can be anything except the empty string.
		TimeKey:        "T",
		LevelKey:       "L",
		NameKey:        "N",
		CallerKey:      "C",
		MessageKey:     "M",
		StacktraceKey:  "S",
		LineEnding:     zapcore.DefaultLineEnding,
		EncodeLevel:    zapcore.CapitalColorLevelEncoder,
		EncodeTime:     TimeEncoder,
		EncodeDuration: zapcore.StringDurationEncoder,
		EncodeCaller:   zapcore.ShortCallerEncoder,
	}
}

// TimeEncoder as name says
func TimeEncoder(t time.Time, enc zapcore.PrimitiveArrayEncoder) {
	enc.AppendString(t.Format("2006-01-02 15:04:05.000"))
}

func setLogger(debug bool, log string) {

	level := zap.InfoLevel
	encoder := NewEncoderConfig()
	if debug {
		level = zap.DebugLevel
		encoder = zap.NewDevelopmentEncoderConfig()
	}

	core := zapcore.NewCore(
		zapcore.NewConsoleEncoder(encoder),
		zapcore.AddSync(os.Stderr),
		level,
	)

	if log != "" {
		w := zapcore.AddSync(&lumberjack.Logger{
			Filename:   log,
			MaxSize:    500, // megabytes
			MaxBackups: 3,
			MaxAge:     28, // days
		})

		core = zapcore.NewCore(
			zapcore.NewJSONEncoder(encoder),
			zapcore.AddSync(w),
			level,
		)
	}

	logger = zap.New(core, zap.AddCaller())
	defer logger.Sync()
	sugar = logger.Sugar()
}
